# Load package
library(edgeR)
library(dplyr)
library(magrittr)

# An useful error if the argument is missing
if (is.null(snakemake@input[["raw_counts"]]) | is.null(snakemake@input[["raw_metadata"]])){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read path of output files
raw_counts_path = snakemake@input[["raw_counts"]]
raw_metadata_path = snakemake@input[["raw_metadata"]]

counts_processed_path = snakemake@output[["counts_processed"]]
metadata_processed_path = snakemake@output[["metadata_processed"]]
genemetadata_path = snakemake@output[["genemetadata"]]
target_ct_file_path = snakemake@output[["target_ct_file"]]

### ----------------------------------------------------------------- ###
##  ---------------------------- Check Data -------------------------  ##
### ----------------------------------------------------------------- ###


counts = read.table(raw_counts_path)
rownames(counts) = toupper(rownames(counts))
metadata = read.table(raw_metadata_path, header = TRUE)
rownames(metadata) = metadata$cell_ID

if(!any("Celltype" == colnames(metadata)) | !any("cell_ID" == colnames(metadata))) 
{
    message("column names of metadata do not contain cell_ID and/or Celltype") ; break
}

# change cell names to remove point
colnames(counts) = gsub("[.-]","_" , colnames(counts))
metadata$cell_ID = gsub("[.-]","_" , metadata$cell_ID)

metadata$Celltype = gsub(" ","." , metadata$Celltype)

# Check if the metadata rows correspond to column names
stopifnot(colnames(counts) == metadata$cell_ID)

target_ct = table(metadata$Celltype) %>% sort(decreasing = T) %>% names %>% extract(1:2)
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
rownames(means_perCT) %<>% toupper

########################################################
# remove genes with mean == 0 in celltypes to simulate #
########################################################

gene_index = which(means_perCT[,target_ct[1]] == 0 | means_perCT[,target_ct[2]] == 0)
if(length(gene_index) > 0) {
  means_perCT = means_perCT[-gene_index,]
  counts = counts[-gene_index,] # remove in counts
  dge = dge[-gene_index,]
}

stopifnot(rownames(means_perCT) == rownames(counts))

######################################
# draw genes from mean-variance plot #
######################################

meanvar_relationship = plotMeanVar(dge, nbins = 20)
df = data.frame(means = meanvar_relationship$bin.means %>% unlist %>% log10, vars = meanvar_relationship$bin.vars %>% unlist %>% log10,
                gene_names = meanvar_relationship$bin.vars %>% unlist %>% names)


# select 10 genes from each bin to be add signal to
set.seed(3)
vec_genes_toAdd_signal = vector()
for(bin in 1:length(meanvar_relationship$bin.means))
{
  n = meanvar_relationship$bin.means[bin] %>% unlist
  names = sample(names(n), size = 10, replace = F)
  vec_genes_toAdd_signal = append(vec_genes_toAdd_signal, names)
}

#############################
# Update gene metadata file #
#############################

genemetadata = list( disp = data.frame(gene = rownames(counts) , edgeR_dispersion = dge$tagwise.dispersion) ,
                      mean = means_perCT %>% as.data.frame )
genemetadata$mean$gene_names = genemetadata$mean %>% rownames()
genemetadata$mean = merge(genemetadata$mean, df, by = "gene_names") # add mean and variance info
genemetadata$mean %<>% mutate(.,gene_to_use = genemetadata$mean$gene_names %in% vec_genes_toAdd_signal)

rownames(genemetadata$disp) = genemetadata$disp$gene %>% toupper

write.table(counts, counts_processed_path , sep = "\t")
write.table(metadata, metadata_processed_path , sep = "\t")
saveRDS(genemetadata, genemetadata_path)
write.table(target_ct, target_ct_file_path , sep = "\t", row.names = F, col.names = F) # sample 2 celtypes with highest amount of cells

sessionInfo()
