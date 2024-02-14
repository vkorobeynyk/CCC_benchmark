# Load package
library(liana)
library(Seurat)
library(dplyr)
library(scuttle)


# An useful error if the argument is missing
if (is.null(snakemake@input[["sc_inflated_counts"]]) | is.null(snakemake@input[["sc_metadata"]]) | is.null(snakemake@input[["simulated_interactions"]]) | is.null(snakemake@output[["significant_interactions"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read the argument
path_sc_inflated_counts <- snakemake@input[["sc_inflated_counts"]]
path_sc_metadata <- snakemake@input[["sc_metadata"]]
path_significant_interactions <- snakemake@output[["significant_interactions"]]
path_simulated_interactions <- snakemake@input[["simulated_interactions"]]

#############
# load data #
#############

raw_counts = read.csv(path_sc_inflated_counts,sep="\t") %>% as.matrix
metadata = read.csv(path_sc_metadata,sep="\t")
simulated_interactions = readRDS(path_simulated_interactions)

#################
# Preprocessing #
#################

SO = CreateSeuratObject(raw_counts, meta.data = metadata)
SO = NormalizeData(SO)
Idents(SO) = metadata$Celltype


##############
# Run method #
##############

# Run liana method
method_out = liana_wrap(SO, method = "sca", resource = "OmniPath", min_cells = 0, return_all = T)

######################
# Create output list #
######################

significant_interactions = vector(mode = 'list', length = length(simulated_interactions))
names(significant_interactions) = names(simulated_interactions)

# add a CT_CT column 
method_out$source_target = paste0(method_out$source , "_" , method_out$target)

# Select only significant interactions
method_out = method_out[method_out$LRscore > 0.5, c("ligand", "receptor", "LRscore" , "source_target")]

# iterate over all CT_CT combinations and append the significant interactions to the list
for(CT_CT in names(significant_interactions))
{
    significant_interactions[[CT_CT]] = filter(method_out ,source_target == CT_CT ) %>% select(. , c("ligand","receptor","LRscore"))
}


# save file
saveRDS(significant_interactions , path_significant_interactions)
