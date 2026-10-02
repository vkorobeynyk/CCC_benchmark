suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(glue)
  library(yaml)
})

# =============================================================================
# 1. Environment and Parameter Setup
# =============================================================================
if (is.null(snakemake@input[["significant_interactions"]]) | 
    is.null(snakemake@input[["genemetadata"]]) |
    is.null(snakemake@output[["final_scores"]])) {
  stop("Missing Snakemake arguments.")
}

path_final_scores = snakemake@output[["final_scores"]]
genemetadata_paths = snakemake@input[["genemetadata"]]
config = yaml::read_yaml("config.yaml")

methods = config$methods %>% unlist()
datasets = config$datasets %>% unlist()
FC_values = config$semiSimulation$FC %>% unlist() %>% as.double()
PCE_values = config$semiSimulation$PCE %>% unlist() %>% as.double()

#' Compute both detection rate and rank quality
#'
#' @param df           Method's full output table for one FC/PCE run (pre-significance-filter)
#' @param target_genes Simulated (ground-truth) ligand_receptor IDs
#'
#' @return List with:
#'   detection_rate         - fraction of target_genes actually present in df
#'   normalizedRank_conditional - mean normalized rank, computed ONLY among
#'                             target genes that were detected (excludes NA)
#'   normalizedRank_penalized   - the original metric: missing genes forced to
#'                             worst-case rank (N)
calculate_normalizedRank = function(df, target_genes) {
  N = nrow(df)
  ranks = match(target_genes, df$ligand_receptor)
  
  is_detected = !is.na(ranks)
  detection_rate = mean(is_detected)
  
  # Conditional rank: only among genes the method actually reported
  norm_ranks_conditional = (ranks[is_detected] - 1) / (N - 1)
  mean_rank_conditional = if (length(norm_ranks_conditional) > 0) {
    mean(norm_ranks_conditional)
  } else {
    NA  # method detected NONE of the target genes -- rank quality is undefined, not 1
  }
  
  # Penalized version: undetected -> worst rank
  ranks_penalized = ranks
  ranks_penalized[is.na(ranks_penalized)] = N
  mean_rank_penalized = mean((ranks_penalized - 1) / (N - 1))
  
  list(
    detection_rate               = detection_rate,
    normalizedRank_conditional   = mean_rank_conditional,
    normalizedRank_penalized     = mean_rank_penalized
  )
}

# =============================================================================
# 2. Main Processing Loop
# =============================================================================
statistics_lst = list()

for (dataset in datasets) {
  
  # Load ground truth (simulated interactions) directly from genemetadata --
  # these are the same for the entire dataset (computed once in
  # processing_dataset.R), so no need to read a separate per-FC/PCE file
  genemetadata_path = genemetadata_paths[grepl(paste0("^data/processed/", dataset, "/"), genemetadata_paths)]
  genemetadata = readRDS(genemetadata_path)
  
  simulated_interactions_lst = genemetadata$simulated_interactions_lst
  
  for (method in methods) {
    
    # Iterate over parameter combinations using a grid
    params_grid = expand.grid(FC = FC_values, PCE = PCE_values)
    
    for (i in 1:nrow(params_grid)) {
      FC = params_grid[i, "FC"]
      PCE = params_grid[i, "PCE"]
      naming = glue("significant_interactions_FC_{FC}_PercCellsExpressing_{PCE}.csv")
      file_path = file.path("output", dataset, method, naming)
      
      # Block-style if for file check
      if (!file.exists(file_path)) {
        break
      }
      
      # Load method output
      returned_interactions_df = read.table(file_path, header = TRUE, sep = "\t")
      
      # Extract interactions flagged as significant
      significant_detected = returned_interactions_df %>% 
        dplyr::filter(significant == "True") %>% 
        pull(ligand_receptor)
      
      # ---------------------------------------------------------------------
      # Metric Calculations
      # ---------------------------------------------------------------------
      
      # Percentage of significant interactions among simulated
      perc_significantInteractions_among_simulated = mean(simulated_interactions_lst$ligand_receptor %in% significant_detected) * 100
      
      tp = sum(significant_detected %in% simulated_interactions_lst$ligand_receptor)
      fp = length(significant_detected) - tp
      fn = length(simulated_interactions_lst$ligand_receptor) - tp
      
      if (tp + fp > 0) {
        precision = tp / (tp + fp)
      } else {
        precision = 0
      }
      
      if (tp + fn > 0) {
        recall = tp / (tp + fn)
      } else {
        recall = 0
      }
      
      if (precision + recall > 0) {
        f1 = 2 * (precision * recall) / (precision + recall)
      } else {
        f1 = 0
      }
      
      norm_rank = calculate_normalizedRank(returned_interactions_df, simulated_interactions_lst$ligand_receptor)
      print(which(simulated_interactions_lst$ligand_receptor %in% returned_interactions_df$ligand_receptor))
      # Store results directly in the nested list structure
      statistics_lst[[dataset]][[method]][[naming]] = data.frame(
        dataset = dataset,
        method = method,
        FC = FC,
        PCE = PCE,
        TP = tp,
        FP = fp,
        FN = fn,
        precision = precision %>% round(3),
        recall = recall %>% round(3),
        f1score = f1 %>% round(3),
        perc_significantInteractions_among_simulated = perc_significantInteractions_among_simulated,
        detection_rate = norm_rank$detection_rate %>% round(3),
        normalizedRank_conditional = norm_rank$normalizedRank_conditional %>% round(3),
        normalizedRank_penalized = norm_rank$normalizedRank_penalized %>% round(3)
      )
    }
  }
}

# =============================================================================
# 3. Save Aggregated Scores
# =============================================================================
saveRDS(statistics_lst, path_final_scores)