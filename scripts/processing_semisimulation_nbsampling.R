suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(edgeR)
})
source("scripts/helper_functions.R")

# =============================================================================
# 1. Environment and Parameter Setup
# =============================================================================
# Inputs
counts_processed_path = snakemake@input[["processed_counts"]]
metadata_processed_path = snakemake@input[["metadata_processed"]]
genemetadata_path = snakemake@input[["genemetadata"]]

# Outputs
path_sc_inflated_counts = snakemake@output[["inflated_counts"]]
path_sc_metadata = snakemake@output[["metadata_processed"]]
path_PCE_perGene = snakemake@output[["PCE_perGene"]]
path_simulated_interactions = snakemake@output[["simulated_interactions"]]
path_FC_after_simulation = snakemake@output[["FC_after_simulation"]]

# Params/Wildcards
nLR_per_CTCTcomb = snakemake@params[["nLR_per_CTCTcomb"]]
LR_database_path = snakemake@params[["LR_database"]]
FC = as.double(snakemake@wildcards[["FC"]])
PCE = as.integer(snakemake@wildcards[["PCE"]])

# =============================================================================
# 2. Data Loading and Filtering
# =============================================================================
counts = read.table(counts_processed_path, header = TRUE, row.names = 1, sep = "\t")
metadata = read.table(metadata_processed_path, header = TRUE, sep = "\t")
genemetadata = readRDS(genemetadata_path)
means_perCT = genemetadata$mean
rownames(means_perCT) = means_perCT$gene_names
original_offsets = genemetadata$edgeR_offsets ; rownames(original_offsets) = original_offsets$cell_ID

identical(metadata$cell_ID, genemetadata$edgeR_offsets$cell_ID) # precaution

'
counts = read.table("/home/vkorob/Documents/git/CCC_benchmark_afterLLMoverhaul//data/processed/10x/counts_10x_processed.tsv")
metadata = read.table("/home/vkorob/Documents/git/CCC_benchmark_afterLLMoverhaul//data/processed/10x/metadata_10x_processed.tsv", header = TRUE, sep = "\t")
genemetadata = readRDS("/home/vkorob/Documents/git/CCC_benchmark_afterLLMoverhaul//data/processed/10x/genemetadata.RDS")
means_perCT = genemetadata$mean
LR_database_path = "data/LR_database.tsv"
'

# Ensure cell alignment
stopifnot(colnames(counts) == metadata$cell_ID)

# =============================================================================
# 3. Load Interaction to be simulated (Ligand-Receptor Database)
# =============================================================================
simulated_interactions_lst = genemetadata$simulated_interactions_lst

# =============================================================================
# 4. Semi-Simulation
# =============================================================================
# Perform NB sampling to inflate counts for selected genes
semi_simulation_out = semi_simulate(
  counts = counts,
  simulated_interactions_lst = simulated_interactions_lst,
  genemetadata = genemetadata,
  metadata = metadata,
  FC = FC,
  PCE = PCE
)

# =============================================================================
# 5. Effective Fold Change Calculation
# =============================================================================
# Re-estimate means to calculate the actual realized Fold Change
dge_inflated = DGEList(counts = semi_simulation_out$counts_inflated, samples = metadata)
FC_after_semisimulation = list()

# Focus only on cells where counts were modified
cells_to_keep = setdiff(colnames(dge_inflated), unlist(semi_simulation_out$not_inflated_cells))
tmp_dge = dge_inflated[, cells_to_keep]
metadata_subset = metadata %>% filter(cell_ID %in% cells_to_keep)
design_matrix = model.matrix(as.formula("~0 + Celltype"), data = metadata_subset)

tmp_dge = estimateDisp(tmp_dge, design = design_matrix)

# Estimate effective mean (mu)
logmeans = edgeR::mglmOneWay(
  tmp_dge$counts, 
  offset = original_offsets[match(colnames(tmp_dge), original_offsets$cell_ID), "edgeR_offset"], # contains original library sizes
  design = design_matrix,
  dispersion = tmp_dge$tagwise.dispersion
)

means_real = exp(logmeans$coefficients) %>% as.data.frame()
colnames(means_real) = gsub("Celltype", "", colnames(design_matrix))
means_real$gene_names = rownames(means_real)

FC_after_semisimulation[["Sender"]] = 
  means_real[semi_simulation_out$baseline_rate_lst$Sender %>% names,"Sender"] /
  semi_simulation_out$baseline_rate_lst$Sender %>% unlist

stopifnot(names(semi_simulation_out$baseline_rate_lst$Sender) == means_real[semi_simulation_out$baseline_rate_lst$Sender %>% names,"gene_names"]) # precaution

FC_after_semisimulation[["Receiver"]] = 
  means_real[semi_simulation_out$baseline_rate_lst$Receiver %>% names,"Receiver"]/
  semi_simulation_out$baseline_rate_lst$Receiver %>% unlist

# =============================================================================
# 6. Save Results
# =============================================================================
write.table(semi_simulation_out$counts_inflated, path_sc_inflated_counts, sep = "\t", quote = FALSE)
write.table(metadata, path_sc_metadata, sep = "\t", quote = FALSE)
saveRDS(semi_simulation_out$PCE_lst, path_PCE_perGene)
saveRDS(FC_after_semisimulation, path_FC_after_simulation)