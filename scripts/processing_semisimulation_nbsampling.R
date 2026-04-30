library(dplyr)
library(stringr)
library(magrittr)
library(edgeR)
source("scripts/helper_functions.R")

# INPUT FILES
counts_processed_path = snakemake@input[["processed_counts"]]
metadata_processed_path = snakemake@input[["metadata_processed"]]
genemetadata_path = snakemake@input[["genemetadata"]]

# OUTPUT FILES
path_sc_inflated_counts = snakemake@output[["inflated_counts"]]
path_sc_metadata = snakemake@output[["metadata_processed"]]
path_PCE_perGene = snakemake@output[["PCE_perGene"]]
path_simulated_interactions = snakemake@output[["simulated_interactions"]]
path_FC_after_simulation = snakemake@output[["FC_after_simulation"]]

##########
# Params #
##########

nLR_per_CTCTcomb = snakemake@params[["nLR_per_CTCTcomb"]]
FC = as.double(snakemake@wildcards[["FC"]])
LR_database_path = snakemake@params[["LR_database"]]
PCE = as.integer(snakemake@wildcards[["PCE"]])

#############
# read data #
#############

counts = read.table(counts_processed_path)
metadata = read.table(metadata_processed_path)
genemetadata = readRDS(genemetadata_path)
means_perCT = genemetadata$mean

'
counts = read.table("/home/vkorob/Documents/git/CCC_benchmark/data/processed/10x/counts_10x_processed.tsv")
metadata = read.table("/home/vkorob/Documents/git/CCC_benchmark/data/processed/10x/metadata_10x_processed.tsv")
genemetadata = readRDS("/home/vkorob/Documents/git/CCC_benchmark/data/processed/10x/genemetadata.RDS")
means_perCT = genemetadata$mean
LR_database_path = "data/LR_database.tsv"
'

# check if cell names of counts and metadata correspond and are in the same order
stopifnot(colnames(counts) == metadata$cell_ID)

# remove genes
message(paste("Shape of count dataframe:" , str_flatten(dim(counts) , " ")))

########################################################
# remove genes with mean == 0 in celltypes to simulate #
########################################################

gene_index = which(means_perCT[,"CT1"] == 0 | means_perCT[,"CT2"] == 0)
if(length(gene_index) > 0) {
  means_perCT = means_perCT[-gene_index,]
}

####################
# Load LR_database #
####################

# Load database
LR_database = read.table(LR_database_path, header = T)

# iterate over every row and filter genes that are not in count data and have 0 mean per CT
n = LR_database %>%
  apply(., 1, function(row) {
    str_split(row,"_") %>%
      lapply(., function(x) {(x %in% rownames(counts)) & (x %in% means_perCT$gene_names)}) %>%
      unlist %>%
      all
  })
LR_database = LR_database[n,]

message(paste("After filtering LR database," , nrow(LR_database) , "LR pairs show expression in at least 10 cells"))

###########################################
# Select genes to artificially add counts #
###########################################
simulated_interactions_lst = list()

# Pre-sample LR pairs to be used in the semi simulation
# Also pre-sample the subunit genes

# select LR pair to inflate expression according to index indexLR_toSample
if(nrow(LR_database) < nLR_per_CTCTcomb) {stop("Amount of LR pairs to sample is higher than available in dataset")}
LR_sample = LR_database[1:nLR_per_CTCTcomb,]
simulated_interactions_lst[["ligand"]] = LR_sample$ligand %>% str_split(., "_") %>% unlist %>% unique
simulated_interactions_lst[["receptor"]] = LR_sample$receptor %>% str_split(., "_") %>% unlist %>% unique
simulated_interactions_lst[["ligand_receptor"]] = LR_sample$ligand_receptor

# Semi simulation
# It may happen that a subunit of a gene has a mean parameter that is below the 1Q threshold I use. Now I use the original means for the subunits. IN theory i would have 
# to check the parameters of every subunit and they are not in line, I would remove them
semi_simulation_out = semi_simulate(counts = counts, simulated_interactions_lst = simulated_interactions_lst , genemetadata = genemetadata, 
                                    metadata = metadata , FC = FC, pce = PCE)

#################################
# calculate FC after simulation #
#################################
dge = DGEList(counts = semi_simulation_out$counts_inflated, samples = metadata)


# estimate the mean only for cells that had their expression increased
# iterate over 2 CT and subset dge to only contain the inflated cells for that CT

FC_after_semisimulation = list()
for(CT in names(semi_simulation_out$not_inflated_cells))
{
  # select the proper set of genes (either ligand for sender Ct or receiver for receiving CT)
  if(CT == "CT1") {simulated_genes = simulated_interactions_lst[["ligand"]]} else if(CT == "CT2") {simulated_genes = simulated_interactions_lst[["receptor"]]}
  
  cells_to_keep = setdiff(colnames(dge), semi_simulation_out$not_inflated_cells)
  tmp_dge = dge[,cells_to_keep]
  
  # update model matrix
  metadata2 = metadata %>% filter(cell_ID %in% cells_to_keep) # filter metadata
  
  stopifnot(metadata2$cell_ID == colnames(tmp_dge)) # just for security
  
  mm= model.matrix(as.formula("~0 + Celltype") , metadata2)
  
  # Estimate disp
  tmp_dge = estimateDisp(tmp_dge , design = mm)
  tmp_dge = edgeR::calcNormFactors(tmp_dge)
  
  # estimating mu
  centered.off = edgeR::getOffset(tmp_dge)  
  logmeans = edgeR::mglmOneWay(tmp_dge$counts, offset = 0, design = mm,
                                dispersion = tmp_dge$tagwise.dispersion) 
  
  means_perCT_semisimulation = exp(logmeans$coefficients) %>% as.data.frame
  colnames(means_perCT_semisimulation) = colnames(mm)
  colnames(means_perCT_semisimulation) = gsub("Celltype","",colnames(means_perCT_semisimulation))
  means_perCT_semisimulation$gene_names = rownames(means_perCT_semisimulation)
  
  FC_after_semisimulation[[CT]] = means_perCT_semisimulation %>% arrange(., gene_names) %>% filter(.,gene_names %in% simulated_genes) %>% .[CT] / 
    means_perCT %>% arrange(., gene_names) %>% filter(.,gene_names %in% simulated_genes) %>% .[CT]
}


################
# save results #
################

write.table(semi_simulation_out$counts_inflated, path_sc_inflated_counts , sep = "\t")
write.table(metadata,path_sc_metadata, sep = "\t")
saveRDS(simulated_interactions_lst,path_simulated_interactions)
saveRDS(semi_simulation_out$PCE_lst,path_PCE_perGene)
saveRDS(FC_after_semisimulation,path_FC_after_simulation)