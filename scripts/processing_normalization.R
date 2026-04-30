suppressMessages({
  library(scater)
  library(SingleCellExperiment)
  library(dplyr)
})

#############
### INPUT ###
#############
inflated_counts_path = snakemake@input[["inflated_counts"]] # Processed and inflated counts counts

##############
### OUTPUT ###
##############
normalized_counts_path = snakemake@output[["inflated_normalized_counts"]] # Processed, inflated and normalized counts

inflated_counts = read.table(inflated_counts_path)

# Normalize counts using scater package
normalized_inflated_counts = SingleCellExperiment(list(counts = inflated_counts)) %>% logNormCounts %>% assay(., "logcounts")

write.table(normalized_inflated_counts , normalized_counts_path, sep = "\t")