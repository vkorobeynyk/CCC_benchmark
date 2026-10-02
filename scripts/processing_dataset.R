suppressMessages({
  library(edgeR)
  library(dplyr)
  library(magrittr)
  library(stringr)
})
source("scripts/helper_functions.R")

# =============================================================================
# 1. Environment and Path Setup
# =============================================================================
if (is.null(snakemake@input[["raw_counts"]]) | is.null(snakemake@input[["raw_metadata"]]) | is.null(snakemake@params[["LR_database"]]) | is.null(snakemake@params[["nLR_per_CTCTcomb"]])) {
  stop("Input paths for some variables are missing.")
}

# Input paths
raw_counts_path = snakemake@input[["raw_counts"]]
raw_metadata_path = snakemake@input[["raw_metadata"]]

#raw_counts_path = "data/10x/raw_counts.tsv"
#raw_metadata_path = "data/10x/raw_metadata.tsv"

# Output paths
counts_processed_path = snakemake@output[["processed_counts"]]
metadata_processed_path = snakemake@output[["metadata_processed"]]
genemetadata_path = snakemake@output[["genemetadata"]]

# Params
LR_database_path = snakemake@params[["LR_database"]]
nLR_per_CTCTcomb = snakemake@params[["nLR_per_CTCTcomb"]]

# =============================================================================
# 2. Data Loading and Nomenclature Standardization
# =============================================================================
counts = read.table(raw_counts_path, header = TRUE, row.names = 1)
metadata = read.table(raw_metadata_path, header = TRUE)
LR_database = read.table(LR_database_path, header = TRUE, sep = " ")

# Standardize gene names to uppercase
rownames(counts) = toupper(rownames(counts))

# Ensure metadata is indexed by cell_ID
rownames(metadata) = metadata$cell_ID

# Validation check for required columns
if (!any("Celltype" == colnames(metadata)) | !any("cell_ID" == colnames(metadata))) {
  stop("Metadata must contain 'cell_ID' and 'Celltype' columns.")
}

# Remove characters that can cause issues in downstream naming (dots and hyphens)
colnames(counts) = gsub("[.-]", "_", colnames(counts))
metadata$cell_ID = gsub("[.-]", "_", metadata$cell_ID)
rownames(metadata) = metadata$cell_ID
metadata$Celltype = gsub(" ", ".", metadata$Celltype)

# Verify that counts and metadata align
stopifnot(colnames(counts) == metadata$cell_ID)

# =============================================================================
# 3. Filtering
# =============================================================================
# Filter genes expressed in at least 10 cells and with a total count > 200
counts = counts[rowSums(counts != 0) > 10, ]
counts = counts[rowSums(counts) > 200, ]

# Update metadata to match the filtered cell set if necessary
metadata = metadata[colnames(counts), ]

# =============================================================================
# 4. Statistical Estimation (edgeR)
# =============================================================================
dge = DGEList(counts = counts, samples = metadata)
design_matrix = model.matrix(as.formula("~0 + Celltype"), data = metadata)

# Estimate dispersion and normalization factors
dge = estimateDisp(dge, design = design_matrix)
dge = edgeR::calcNormFactors(dge)

# Estimate mean expression (mu) per cell type
offset = edgeR::getOffset(dge)
logmeans = edgeR::mglmOneWay(
  dge$counts, 
  offset = offset, 
  design = design_matrix,
  dispersion = dge$tagwise.dispersion
)

means_perCT = exp(logmeans$coefficients) %>% as.data.frame
colnames(means_perCT) = gsub("Celltype", "", colnames(design_matrix))
rownames(means_perCT) = toupper(rownames(means_perCT))
means_perCT = means_perCT %>% mutate(gene_names = rownames(.))

# =============================================================================
# 5. Determine Rate Threshold for LR Pair Selection
# =============================================================================
# This step is important as we dont want to select for simulation LR genes
# that are very sparsely expressed (1%) of cells as their estimated means are too small

# Different sequencing technologies vary widely in sparsity, so a fixed rate
# threshold doesn't yield a comparable number of usable LR pairs across
# datasets. Here we search for the strictest threshold that still yields at
# least nLR_per_CTCTcomb valid pairs for this specific dataset, and store it
# in genemetadata for use downstream in the semi-simulation script.

threshold_result = find_threshold_for_target_pairs(
  LR_database = LR_database,
  means_perCT = means_perCT,
  offset      = offset,
  metadata    = metadata,
  counts      = counts,
  target_n    = nLR_per_CTCTcomb
)

# Specify the threshold that yields at least 30 pairs
LR_rate_threshold = threshold_result$threshold

# Select top N interactions and extract unique genes for ligands and receptors
LR_database = LR_database[threshold_result$valid_rows, ]

if (nrow(LR_database) < nLR_per_CTCTcomb) {
  stop("Insufficient valid LR pairs after rate filtering.")
}
LR_sample = LR_database[1:nLR_per_CTCTcomb, ]
simulated_interactions_lst = list(
  ligand = LR_sample$ligand %>% str_split("_") %>% unlist() %>% unique(),
  receptor = LR_sample$receptor %>% str_split("_") %>% unlist() %>% unique(),
  ligand_receptor = LR_sample$ligand_receptor
)

# =============================================================================
# 6. Mean-Variance Relationship and Random Sampling
# =============================================================================
# Generate mean-variance relationship for diagnostic visualization
meanvar_relationship = plotMeanVar(dge, nbins = 20)
df_meanvar = data.frame(
  means = unlist(meanvar_relationship$bin.means) %>% log10(),
  vars = unlist(meanvar_relationship$bin.vars) %>% log10(),
  gene_names = names(unlist(meanvar_relationship$bin.vars))
)

# Sample 10 random genes from each bin to serve as a baseline for the simulation
set.seed(3)
vec_random_genes = vector()
for (bin in 1:length(meanvar_relationship$bin.means)) {
  bin_genes = names(unlist(meanvar_relationship$bin.means[bin]))
  
  if (length(bin_genes) >= 10) {
    sampled_names = sample(bin_genes, size = 10, replace = FALSE)
  } else {
    sampled_names = bin_genes
  }
  vec_random_genes = append(vec_random_genes, sampled_names)
}

# =============================================================================
# 7. Construct and Save Metadata
# =============================================================================
genemetadata = list(
  disp = data.frame(
    gene = rownames(dge), 
    edgeR_dispersion = dge$tagwise.dispersion
  ),
  mean = means_perCT,
  edgeR_offsets = data.frame(cell_ID = colnames(dge), edgeR_offset = offset),
  LR_rate_threshold = LR_rate_threshold,
  simulated_interactions_lst = simulated_interactions_lst
)

# make sure cell order is same
stopifnot(genemetadata$edgeR_offsets$cellnames == metadata$cell_ID)
stopifnot(colnames(counts) == metadata$cell_ID)
# Add auxiliary info to mean metadata
genemetadata$mean = merge(genemetadata$mean, df_meanvar, by = "gene_names")
genemetadata$mean$randomly_sampled = genemetadata$mean$gene_names %in% vec_random_genes

# Sort gene names
genemetadata$disp = genemetadata$disp %>% arrange(gene)
genemetadata$mean = genemetadata$mean %>% arrange(gene_names)

# Final export
write.table(counts, counts_processed_path, sep = "\t", quote = FALSE)
write.table(metadata, metadata_processed_path, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(genemetadata, genemetadata_path)