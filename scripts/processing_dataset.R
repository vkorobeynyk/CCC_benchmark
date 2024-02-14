# Load package
library(edgeR)
library(dplyr)

# An useful error if the argument is missing
if (is.null(snakemake@input[["raw_counts"]]) | is.null(snakemake@input[["raw_metadata"]])){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read path of output files
raw_counts_path = snakemake@input[["raw_counts"]]
raw_metadata_path = snakemake@input[["raw_metadata"]]

counts_processed_path = snakemake@output[["counts_processed"]]
metadata_processed_path = snakemake@output[["metadata_processed"]]
gene_metadata_path = snakemake@output[["gene_metadata"]]
target_ct_file_path = snakemake@output[["target_ct_file"]]

### ----------------------------------------------------------------- ###
##  ---------------------------- Check Data -------------------------  ##
### ----------------------------------------------------------------- ###


counts = read.table(raw_counts_path)
metadata = read.table(raw_metadata_path, header = TRUE)

if(!any("Celltype" == colnames(metadata)) | !any("cell_ID" == colnames(metadata))) 
{
    message("column names of metadata do not contain cell_ID and/or Celltype") ; break
}

# change cell names to remove point
colnames(counts) = gsub("[.]","_" , colnames(counts))
metadata$cell_ID = gsub("[-]","_" , metadata$cell_ID)

metadata$Celltype = gsub(" ","." , metadata$Celltype)


# Check if the metadata rows correspond to column names
stopifnot(colnames(counts) == metadata$cell_ID)

#############
# Filtering #
#############

# Filter genes based on minimal amount of cells expressing the gene
counts = counts[rowSums(counts != 0) > 10,]

# Filter genes based on minimal amount of counts it should express
counts = counts[rowSums(counts) > 200,]


message(paste("Filtering - gene has to be expressed in at least 10 cells and have 200 counts"))

######################################
# estimate mean and disp using edgeR #
######################################

dge = DGEList(counts = counts, samples = metadata)
mm= model.matrix(as.formula("~0 + Celltype") , metadata)

# Estimate disp
dge <- estimateDisp(dge , design = mm)
dge <- edgeR::calcNormFactors(dge)

# estimating mu
centered.off <- edgeR::getOffset(dge)  
centered.off <- centered.off - mean(centered.off) 
logmeans <- edgeR::mglmOneWay(dge$counts, offset = centered.off, design = mm,
                                dispersion = dge$tagwise.dispersion) 

means_perCT = exp(logmeans$coefficients)
colnames(means_perCT) = colnames(mm)
colnames(means_perCT) = gsub("Celltype","",colnames(means_perCT))

#############################
# Update gene metadata file #
#############################
gene_metadata = list( disp = data.frame(gene = rownames(counts) , edgeR_dispersion = dge$tagwise.dispersion) ,
                      mean = means_perCT)

write.table(counts, counts_processed_path , sep = "\t")
write.table(metadata, metadata_processed_path , sep = "\t")
saveRDS(gene_metadata, gene_metadata_path)
set.seed(1)
write.table(sample(table(metadata$Celltype),4) %>% names, target_ct_file_path , sep = "\t", row.names = F, col.names = F) # sample 4 celltypes

sessionInfo()
