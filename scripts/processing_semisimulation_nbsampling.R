library(dplyr)
library(stringr)
library(liana)
library(magrittr)
library(edgeR)
source("scripts/helper_functions.R")

# An useful error if the argument is missing
if (is.null(snakemake@input[["counts_processed"]]) | is.null(snakemake@input[["metadata_processed"]]) | is.null(snakemake@input[["gene_metadata"]]) | is.null(snakemake@input[["target_ct_file"]]) | 
    is.null(snakemake@params[["nLR_per_CTCTcomb"]]) | is.null(snakemake@wildcards[["FC"]]) | is.null(snakemake@wildcards[["perc_cells_expressing"]]) | 
    is.null(snakemake@output[["sc_inflated_counts"]]) | is.null(snakemake@output[["simulated_interactions"]])  | is.null(snakemake@output[["sc_metadata"]])  | is.null(snakemake@output[["perc_cells_expressing_perGene"]])){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}
# OUTPUT FILES
path_sc_inflated_counts <- snakemake@output[["sc_inflated_counts"]]
path_sc_metadata <- snakemake@output[["sc_metadata"]]
path_perc_cells_expressing_perGene <- snakemake@output[["perc_cells_expressing_perGene"]]
path_simulated_interactions <- snakemake@output[["simulated_interactions"]]
path_FC_after_simulation <- snakemake@output[["FC_after_simulation"]]


# INPUT FILES
counts_processed_path <- snakemake@input[["counts_processed"]]
metadata_processed_path <- snakemake@input[["metadata_processed"]]
gene_metadata_path <- snakemake@input[["gene_metadata"]]
target_ct_file_path <- snakemake@input[["target_ct_file"]]

##################
# Set parameters #
##################

nLR_per_CTCTcomb <- snakemake@params[["nLR_per_CTCTcomb"]]
FC <- as.double(snakemake@wildcards[["FC"]])
perc_cells_expressing <- as.integer(snakemake@wildcards[["perc_cells_expressing"]])


# output files
sc_inflated_counts_path <- snakemake@output[["sc_inflated_counts"]]
simulated_interactions_path <- snakemake@output[["simulated_interactions"]]
sc_metadata_path <- snakemake@output[["sc_metadata"]]
cells_sampled_perCTCT_path <- snakemake@output[["cells_sampled_perCTCT"]]

#############
# read data #
#############

counts = read.table(counts_processed_path)
rownames(counts) = toupper(rownames(counts))
metadata = read.table(metadata_processed_path)
rownames(metadata) = metadata$cell_ID
genemetadata = readRDS(gene_metadata_path)
rownames(genemetadata$disp) = genemetadata$disp$gene %>% toupper
means_perCT = genemetadata$mean
rownames(means_perCT) = rownames(means_perCT) %>% toupper
target_ct = read.table(target_ct_file_path) %>% unlist %>% as.character

#counts = read.table("/home/vkorob/Documents/git/CCC_benchmark/data/processed/VASAseq/counts_VASAseq_processed.tsv")
#rownames(counts) = toupper(rownames(counts))
#metadata = read.table("/home/vkorob/Documents/git/CCC_benchmark/data/processed/VASAseq/metadata_VASAseq_processed.tsv")
#rownames(metadata) = metadata$cell_ID
#genemetadata = readRDS("/home/vkorob/Documents/git/CCC_benchmark/data/processed/VASAseq/gene_metadata.tsv")
#rownames(genemetadata$disp) = genemetadata$disp$gene %>% toupper
#means_perCT = genemetadata$mean
#rownames(means_perCT) = rownames(means_perCT) %>% toupper
#target_ct = read.table("/home/vkorob/Documents/git/CCC_benchmark/data/processed/VASAseq//target_ct_file.tsv") %>% unlist %>% as.character


# check if cell names of counts and metadata correspond and are in the same order
stopifnot(colnames(counts) == metadata$cell_ID)

# Subset counts and metadata to contain only the selected CT
#n = which(metadata$Celltype %in% target_ct)
#counts = counts[,n]
#metadata = metadata[n,]

# remove genes
message(paste("Target celltypes:" , str_flatten(target_ct, " ")))
message(paste("Shape of count dataframe:" , str_flatten(dim(counts) , " ")))


#############
# Load LRdb #
#############

# Load OmniPath database
LRdb = select_resource(c('OmniPath'))[[1]]
colnames(LRdb)[1:2] = c("ligand", "receptor")

LRdb = LRdb[(LRdb$ligand %in% rownames(counts)),]
LRdb = LRdb[(LRdb$receptor %in% rownames(counts)),]
# filter LR because there are duplicated pairs (only 1)
LRdb$L_R = str_c(LRdb$ligand,"_",LRdb$receptor)
LRdb = LRdb[!duplicated(LRdb$L_R),]

message(paste("After filtering LR database," , nrow(LRdb) , "LR pairs show expression in at least 10 cells"))

###########################
# Inflate gene expression #
###########################
set.seed(3)
#combination_CT = expand.grid(target_ct,target_ct)
#combination_CT = paste0(combination_CT$Var1, "_", combination_CT$Var2)
combination_CT = str_flatten(target_ct,"_")
simulated_interactions_lst = list()

tmp_LRdb = LRdb

# Pre-sample LR pairs to be used in the semi simulation
# Also pre-sample the subunit genes
for(comb_CT in combination_CT)
{
  # select randomly LR pair to inflate expression
  index_toSample = sample(seq(1,nrow(tmp_LRdb)) , size = nLR_per_CTCTcomb, replace = F)
  LR_sample = tmp_LRdb[index_toSample,]
  
  # Add genes belonging to a complex to be also inflated (based on CellPhoneDB and CellChatDB DB)
  # CellChat
  LRdb_cellchat = select_resource(c('CellChatDB'))[[1]]
  LRdb_cellchat = rbind(filter(LRdb_cellchat, grepl("COMPLEX", source)) , filter(LRdb_cellchat, grepl("COMPLEX", target)))
  
  # CellPhoneDB
  LRdb_cpdb = select_resource(c('CellPhoneDB'))[[1]]
  LRdb_cpdb = rbind(filter(LRdb_cpdb, grepl("COMPLEX", source)) , filter(LRdb_cpdb, grepl("COMPLEX", target)))
  
  LRdb_join = rbind(LRdb_cellchat , LRdb_cpdb)
  LRdb_join = LRdb_join[!str_c(LRdb_join$source, LRdb_join$target) %>% duplicated,] # remove duplicated entries
  
  # remove the sampled LR pairs 
  tmp_LRdb = tmp_LRdb[-index_toSample,]
  
  # adding L and R to the LR_list to track inflated genes without removing duplicates
  simulated_interactions_lst[[comb_CT]] = LR_sample$L_R
  
  ############## find subunits for ligands
  for(gene in LR_sample$ligand)
  {
    tmp_df = filter(LRdb_join , grepl(gene, source_genesymbol) & grepl("COMPLEX", source))
    
    # in case this gene has no subunits, skip
    if(nrow(tmp_df) == 0) {next}

    subunits = c(tmp_df$source_genesymbol %>% str_split("_") %>% lapply("[",1) , tmp_df$source_genesymbol %>% str_split("_") %>% lapply("[",2)) %>% unlist %>% unique()
    subunits = subunits[!grepl(gene, subunits)] # remove original gene
    
    simulated_interactions_lst[[comb_CT]]  %<>% append(. , subunits[which(subunits %in% rownames(counts))] %>% str_c(., "_subunit")) # remove empty strings and add subunit . Also here we filter subunits that are not present in count data
  }
  
  ############## find subunits for receptors
  for(gene in LR_sample$receptor)
  {
    tmp_df = filter(LRdb_join , grepl(gene, target_genesymbol) & grepl("COMPLEX", target))
    
    # in case this gene has no subunits, skip
    if(nrow(tmp_df) == 0) {next}
    
    subunits = c(tmp_df$target_genesymbol %>% str_split("_") %>% lapply("[",1) , tmp_df$target_genesymbol %>% str_split("_") %>% lapply("[",2)) %>% unlist %>% unique()
    subunits = subunits[!grepl(gene, subunits)] # remove original gene
    
    simulated_interactions_lst[[comb_CT]]  %<>% append(. , subunits[which(subunits %in% rownames(counts))] %>% str_c("subunit_" , .)) # remove empty strings and add subunit . Also here we filter subunits that are not present in count data
  }
  
  # remove duplicates in subunits info
  simulated_interactions_lst[[comb_CT]] = simulated_interactions_lst[[comb_CT]][!duplicated(simulated_interactions_lst[[comb_CT]])]
}

# Semi simulation
semi_simulation_out = semi_simulate(counts = counts, simulated_interactions_lst = simulated_interactions_lst , genemetadata = genemetadata, 
                                    metadata = metadata , combination_CT = combination_CT, FC = FC, pce = perc_cells_expressing)

#################################
# calculate FC after simulation #
#################################

dge = DGEList(counts = semi_simulation_out$counts_inflated, samples = metadata)
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

# get simulated genes
names1 = str_split(simulated_interactions_lst[[1]], "_")  %>% lapply("[[",1) %>% unlist %>% setdiff(.,"subunit") %>% unique
names2 = str_split(simulated_interactions_lst[[1]], "_")  %>% lapply("[[",2) %>% unlist %>% setdiff(.,"subunit") %>% unique

FC_after_semisimulation = list(FC_ct1 = means_perCT[names1,target_ct[1]] / genemetadata$mean[names1,target_ct[1]], 
                                FC_ct2 = means_perCT[names2,target_ct[2]] / genemetadata$mean[names2,target_ct[2]])

################
# save results #
################

write.table(semi_simulation_out$counts_inflated, sc_inflated_counts_path , sep = "\t")
write.table(metadata,sc_metadata_path, sep = "\t")
saveRDS(simulated_interactions_lst,path_simulated_interactions)
saveRDS(semi_simulation_out$perc_cells_expressing_lst,path_perc_cells_expressing_perGene)
saveRDS(FC_after_semisimulation,path_FC_after_simulation)

sessionInfo()
