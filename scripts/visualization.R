suppressMessages({
  library(dplyr)
  library(optparse)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  library(purrr)
  library(magrittr)
  library(plyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(patchwork)
  library(reshape2)
  library(edgeR)
  library(yaml)
  source("scripts/helper_functions.R")
})

# =============================================================================
# 1. Command Line Argument Parsing
# =============================================================================
parse_cli_args = function() {
  option_list = list(
    make_option(c("--config.yaml"), type = "character", default = NULL,
                help = "Path to the YAML configuration file", metavar = "character"),
    make_option(c("--path_output_dir"), type = "character", default = NULL,
                help = "Path to the directory containing simulation output", metavar = "character"),
    make_option(c("--path_results_dir"), type = "character", default = NULL,
                help = "Path to where the final PDF results should be saved", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
  if (is.null(opt$config.yaml) | is.null(opt$path_output_dir) | is.null(opt$path_results_dir)) {
    print_help(opt_parser)
    stop("Required arguments are missing. Execution halted.\n", call. = FALSE)
  }
  
  opt
}

opt = parse_cli_args()
path_config_yaml  = opt$config.yaml
path_output_dir   = opt$path_output_dir
path_results_dir  = opt$path_results_dir

if (!dir.exists(path_results_dir)) {
  dir.create(path_results_dir, recursive = TRUE)
}

# =============================================================================
# 2. Data Loading and Parameter Setup
# =============================================================================
metric_results = readRDS(file.path(path_output_dir, "final_scores.RDS"))
config = yaml::read_yaml(path_config_yaml)

methods  = config$methods  %>% unlist()
datasets = config$datasets %>% unlist()
FC  = config$semiSimulation$FC  %>% unlist() %>% as.double()
PCE = config$semiSimulation$PCE %>% unlist() %>% as.integer()

# =============================================================================
# 3. MAIN: Generate diagnostic plots for each dataset
# =============================================================================
# These lists store diagnostic plots and dataframes across all datasets
diagnostic_plots_perCT       = list()
diagnostic_plots_MeanVar     = list()
diagnostic_plots_realFC      = list()
diagnostic_df_realFC         = list()
diagnostic_gene_densityPlots = list()
simulated_interactions_lst   = list()

for (dataset in datasets) {
  
  file_path = file.path(path_output_dir, paste0(dataset, "_semiSimulation_NB"))
  inflated_counts_files        = list.files(file_path, pattern = "sc_inflated_counts")
  PCE_files                    = list.files(file_path, pattern = "PCE_perGene")
  
  genemetadata = readRDS(file.path("data/processed", dataset, "genemetadata.RDS"))
  simulated_interactions_lst[[dataset]] = genemetadata$simulated_interactions_lst
  
  original_counts = read.table(file.path("data", dataset, "raw_counts.tsv"), header = TRUE, row.names = 1)
  rownames(original_counts) = toupper(rownames(original_counts))
  colnames(original_counts) = gsub("[.]", "_", colnames(original_counts))
  
  master_lst = list()
  
  metadata_path = file.path("data/processed", dataset)
  metadata_file = list.files(metadata_path, pattern = "metadata_")
  metadata = read.table(file.path(metadata_path, metadata_file), header = TRUE)
  
  CT_present  = c("Sender", "Receiver")
  params_grid = expand.grid(FC = FC, PCE = PCE)
  
  # ---------------------------------------------------------------------
  # Loop over the FC/PCE parameter grid for this dataset
  # ---------------------------------------------------------------------
  fc_comparison_df_list = list()
  for (i in 1:nrow(params_grid)) {
    FC_val  = params_grid[i, "FC"]
    PCE_val = params_grid[i, "PCE"]
    
    loaded = load_param_combination(file_path, inflated_counts_files, PCE_files,
                                    simulated_interactions_lst[[dataset]], FC_val, PCE_val)
    param_naming    = loaded$param_naming
    inflated_counts = loaded$inflated_counts
    
    simulated_interactions_lst[[dataset]][[param_naming]] = loaded$ligand_receptor
    master_lst[["PCE"]][[param_naming]] = loaded$realized_pce_pct
    
    # --- Track which genes were inflated in each cell type ---
    inflated_genes = extract_inflated_genes(loaded$ligand_receptor)
    master_lst[["Sender"]][[param_naming]]   = list(L = inflated_genes$ligands)
    master_lst[["Receiver"]][[param_naming]] = list(R = inflated_genes$receptors)
    
    if(param_naming %in% c("FC_1_PCE_10", "FC_5_PCE_20"))
    {
      # --- Gene density plots per cell type ---
      for (CT in CT_present) {
        genes_for_CT = if (CT == "Sender") inflated_genes$ligands else inflated_genes$receptors
        diagnostic_gene_densityPlots[[dataset]][[param_naming]][[CT]] = generate_gene_density_plots(
          genes_for_CT, CT, inflated_counts, original_counts, metadata, FC_val, PCE_val, n_genes = 4
        )
      }
    }
      
    # --- Effective FC diagnostics ---
    fc_comparison_df_list[[param_naming]] = extract_effective_fc(file_path, FC_val, PCE_val)
    
    # --- AveLogCPM (kept instead of full inflated matrices, to save memory) ---
    target_cell_idx = metadata$Celltype %in% CT_present
    avelogcpm = inflated_counts[, target_cell_idx] %>% edgeR::aveLogCPM()
    names(avelogcpm) = rownames(inflated_counts)
    master_lst[["counts_aveLogCPM"]][[param_naming]] = avelogcpm
  }
  
  fc_result = plot_FCafter_semisimulation(do.call(rbind, fc_comparison_df_list), dataset)
  diagnostic_plots_realFC[[dataset]] = fc_result$plot
  diagnostic_df_realFC[[dataset]]    = fc_result$summary
  
  # ---------------------------------------------------------------------
  # Dataset-Level Mean-Variance Plot
  # ---------------------------------------------------------------------
  # the same curated ligand/receptor genes are simulated across the whole FC/PCE
  # grid, so any param_naming gives the same genes
  diagnostic_plots_MeanVar[[dataset]] = build_meanvar_plot(
    genemetadata,
    simulated_ligands   = master_lst[["Sender"]][[param_naming]]   %>% unlist() %>% as.character(),
    simulated_receptors = master_lst[["Receiver"]][[param_naming]] %>% unlist() %>% as.character(),
    dataset = dataset
  )
  
  # ---------------------------------------------------------------------
  # Final Diagnostic Plots (AveLogCPM + PCE consistency, combined across grid)
  # ---------------------------------------------------------------------
  original_counts_diagnostic = original_counts[rownames(inflated_counts), ]
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(
    counts     = original_counts_diagnostic,
    master_lst = master_lst,
    FC_param   = c(min(FC), FC[FC <= median(FC)] %>% tail(1), max(FC)),
    PCE_param  = c(min(PCE), PCE[PCE <= median(PCE)] %>% tail(1), max(PCE)),
    dataset    = dataset,
    metadata   = metadata,
    CT_toPlot  = CT_present
  )
  
  
  # --- Export this dataset's PDF ---
  saveRDS(master_lst, file.path(path_results_dir, paste0(dataset, "_master_lst.RDS")))
  
  pdf(file.path(path_results_dir, paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  ggarrange(plotlist = diagnostic_plots_perCT[[dataset]]$avelogcpm_fixedPCE, common.legend = TRUE) %>% print()
  ggarrange(plotlist = diagnostic_plots_perCT[[dataset]]$PCE_fixedPCE) %>% print()
  
  # Effective fold-change diagnostics, one arranged page per theoretical FC
  print(diagnostic_plots_realFC[[dataset]])
  
  # Mean-variance profile
  print(diagnostic_plots_MeanVar[[dataset]])
  
  # Gene density plots for a fixed set of representative FC/PCE benchmark points
  plot_density_keys = c("FC_1_PCE_10", "FC_5_PCE_20")
  
  for (key in plot_density_keys) {
    plots_for_key = diagnostic_gene_densityPlots[[dataset]][[key]][[1]] # plot only ligands
    if (!is.null(plots_for_key)) {
      ggarrange(plotlist = plots_for_key, common.legend = TRUE) %>% print()
    }
  }
  
  dev.off()
}


# =============================================================================
# 4. MAIN: Precision/recall & significance-overlap plots, per dataset
# =============================================================================
statistics_results_lst_recallprecision_plot = compute_recallprecision_stats(metric_results, diagnostic_df_realFC)

significant_interactions_lst = list()

for (dataset in datasets) {
  
  # --- Load per-method significant interactions for this dataset ---
  for (method in methods) {
    significant_interactions_files = file.path(path_output_dir, dataset, method) %>%
      list.files(pattern = "significant_interactions")
    
    for (file in significant_interactions_files) {
      filepath = file.path(path_output_dir, dataset, method, file)
      
      significant_interactions_lst[["all_significant"]][[method]][[file]] =
        read_significant_ligand_receptors(filepath)
      
      significant_interactions_lst[["intersect_significant_simulated"]][[method]][[file]] =
        intersect(significant_interactions_lst[["all_significant"]][[method]][[file]],
                  simulated_interactions_lst[[dataset]]$ligand_receptor)
    }
  }
  
  # --- F1-score heatmap ---
  f1_heatmap_result = build_f1_heatmap(statistics_results_lst_recallprecision_plot[[dataset]])
  heatmap1 = f1_heatmap_result$heatmap
  
  # --- Per-FC/PCE-combination UpSet plots + binary membership tables ---
  grid_upset_data = build_grid_upset_data(significant_interactions_lst, FC, PCE)
  lst_upset_plots_FC  = grid_upset_data$lst_upset_plots_FC
  lst_upset_plots_PCE = grid_upset_data$lst_upset_plots_PCE
  lst_binarydf        = grid_upset_data$lst_binarydf
  
  # --- % of simulated interactions retrieved, heatmap ---
  total_simulated = config$nLR_per_CTCTcomb
  heatmap2 = build_retrieval_pct_heatmap(
    lst_binarydf, total_simulated,
    ref_dims = list(ncol = ncol(f1_heatmap_result$mat), nrow = nrow(f1_heatmap_result$mat))
  )
  
  # --- F1 dot-plot alternative to the heatmap ---
  #heatmap1_alternative = build_f1_dotplot(statistics_results_lst_recallprecision_plot[[dataset]])
  
  # --- UpSet plot: all significant interactions across every FC/PCE combo ---
  all_significant_interactions = build_method_membership_matrix(
    significant_interactions_lst[["all_significant"]],
    universe = unique(unlist(significant_interactions_lst[["all_significant"]]))
  )
  colnames(all_significant_interactions) = names(significant_interactions_lst[["all_significant"]])
  upset1 = build_upset_plot(all_significant_interactions,
                            "All significant interactions retrieved across all combination of parameters")
  
  # --- UpSet plot: simulated interactions retrieved across every FC/PCE combo ---
  intersect_significant_simulated_interactions = build_method_membership_matrix(
    significant_interactions_lst[["intersect_significant_simulated"]],
    universe = simulated_interactions_lst[[dataset]]$ligand_receptor
  )
  colnames(intersect_significant_simulated_interactions) =
    names(significant_interactions_lst[["intersect_significant_simulated"]])
  upset2 = build_upset_plot(intersect_significant_simulated_interactions,
                            "All simulated interactions retrieved across all combination of parameters")
  
  # --- Understand why some simulated ligand-receptor pairs are never found ---
  intersect_significant_simulated_interactions[rowSums(intersect_significant_simulated_interactions) == 0,]
  
  # --- Precision-recall curves, faceted by PCE and by FC ---
  lst_precision_recall_byPCE = build_precision_recall_plots(statistics_results_lst_recallprecision_plot[[dataset]], "PCE")
  lst_precision_recall_byFC  = build_precision_recall_plots(statistics_results_lst_recallprecision_plot[[dataset]], "FC")
  
  # --- Report rank of simulated interactions in the output of each method
  # --- normalizedRank_penalized -> Normalized rank. if LR is not present in output, then it takes the worst rank
  # --- detection_rate -> % of simulated interactions that are in the method output
  p1 = ggplot(statistics_results_lst_recallprecision_plot[[dataset]], aes(x = FC, y = normalizedRank_penalized, color = method, alpha = detection_rate)) +
    geom_point(size = 2) +
    geom_line(aes(group = method)) +
    facet_grid(~PCE) +
    ggtitle("Ranking of simulated LR interactions") +
    theme_light()
  
  # --- F1-score and normalized-rank summary plots ---
  summary_plots = build_f1_and_rank_plots(statistics_results_lst_recallprecision_plot[[dataset]])
  
  # --- Comparison with baseline runs (no signal added) ---
  baseline_plots = build_baseline_plots(metric_results[[dataset]], dataset)
  
  # --- Save all plots for this dataset ---
  pdf(file.path(path_results_dir, paste0(dataset, "_results_plots.pdf")), width = 15, height = 7)
  wrap_plots(lst_precision_recall_byPCE, guides = "collect") %>% print
  wrap_plots(lst_precision_recall_byFC, guides = "collect")  %>% print
  summary_plots$f1_by_FC   %>% print
  summary_plots$f1_by_PCE  %>% print
  summary_plots$rank_by_PCE %>% print
  p1 %>% print
  heatmap1 %>% print
  #heatmap1_alternative %>% print
  heatmap2 %>% print
  
  lst_upset_plots_FC$FC_1$PCE_20  %>% print
  lst_upset_plots_FC$FC_7$PCE_20 %>% print
  
  baseline_plots$fate        %>% print
  baseline_plots$composition %>% print
  
  upset1 %>% print
  upset2 %>% print
  
  dev.off()
}


# =============================================================================
# 5. MAIN: F1-score comparison across datasets
# =============================================================================
rel_f1  = build_relative_f1_matrix(statistics_results_lst_recallprecision_plot)
x       = rel_f1$x
mat_rel = rel_f1$mat_rel

# NOTE: the original script assigns method_cols three times in a row (Set2
# brewer palette, then Okabe-Ito colorblind-safe palette, then viridis) --
# only the LAST assignment (viridis) actually takes effect for the final
# plot. Preserved here exactly, including the two overwritten assignments,
# to match original output; the first two lines are effectively dead code.
method_cols = setNames(brewer.pal(8, "Set2"), unique(x$method))

okabe_ito_cols = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#000000")
method_cols = setNames(okabe_ito_cols[1:length(unique(x$method))], unique(x$method))

library(viridis)
method_cols = setNames(viridis(8, option = "D"), unique(x$method))

ht = build_relative_f1_stacked_heatmap(mat_rel, method_cols)

png(file.path(path_results_dir, "results_acrossdatasets_plots.png"),
    width = 13, height = 10, units = "in", res = 300)
draw(ht, annotation_legend_list = list(
  Legend(
    labels     = names(method_cols),
    title      = "Methods",
    legend_gp  = gpar(fill = method_cols),
    title_gp   = gpar(fontsize = 12), 
    labels_gp  = gpar(fontsize = 12),
    grid_height = unit(5, "mm"),
    grid_width  = unit(5, "mm")
  )
))
dev.off()

# =============================================================================
# 7. Method consensus rank by technology and expression tier
#
# 1. For each dataset x PCE tier, every method gets its mean score over all FC
#    and PCE values in that tier and is ranked within the dataset (1 = best,
#    ties averaged). A method that retrieves nothing gets the last rank; tiers
#    where no method retrieves anything are left unranked (NA).
# 2. Per technology class x tier, ranks are averaged across datasets, and the
#    methods are ranked again on that mean -> consensus rank (1 = best).
#   cell colour/number = consensus rank
#   brackets           = worst rank of the method in any single dataset
#   "–"                = the method retrieved nothing in any dataset of that class
# =============================================================================
rank_metric = "recall"   # "recall" or "f1score" (precision is ~1 after baseline correction)
n_methods   = length(methods)

technology_class = c(SMARTseq2 = "Plate-based", VASAseq = "Plate-based",
                     `10x` = "Droplet-based", `10x_immune_R1` = "Droplet-based",
                     BD_wholeTranscriptome_P1 = "Droplet-based")

pce_tier = function(pce) {
  cut(pce, breaks = c(-Inf, 9, 25, Inf),
      labels = c("Low (4-7%)", "Medium (10-20%)", "High (30-40%)"))
}

rank_df = do.call(rbind.data.frame, statistics_results_lst_recallprecision_plot) %>%
  dplyr::mutate(PCE   = as.numeric(as.character(PCE)),
                FC    = as.numeric(as.character(FC)),
                score = .data[[rank_metric]],
                tier  = pce_tier(PCE),
                technology = technology_class[dataset]) %>%
  # mean score per dataset x tier x method
  dplyr::group_by(technology, dataset, tier, method) %>%
  dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  # rank within dataset x tier
  dplyr::group_by(technology, dataset, tier) %>%
  dplyr::mutate(rank = dplyr::case_when(
    all(score == 0) ~ NA_real_,                                  # nobody retrieves anything: unranked
    score == 0      ~ as.numeric(n_methods),                     # retrieved nothing: last place
    TRUE            ~ rank(-score, ties.method = "average"))) %>%
  dplyr::ungroup()

rank_summary = rank_df %>%
  dplyr::group_by(technology, tier, method) %>%
  dplyr::summarise(mean_rank   = mean(rank, na.rm = TRUE),
                   worst_rank  = suppressWarnings(max(rank, na.rm = TRUE)),
                   n_datasets  = sum(!is.na(rank)),
                   n_retrieved = sum(score > 0),
                   .groups = "drop") %>%
  dplyr::mutate(mean_rank  = ifelse(n_datasets == 0, NA, mean_rank),
                worst_rank = ifelse(n_datasets == 0, NA, worst_rank)) %>%
  # consensus rank within each technology x tier (methods that retrieved nothing are not ranked)
  dplyr::group_by(technology, tier) %>%
  dplyr::mutate(consensus_rank = rank(ifelse(n_retrieved == 0, NA, mean_rank),
                                      ties.method = "min", na.last = "keep")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(column = paste(technology, tier, sep = " | "))

# --- matrices: methods x (technology | tier) ---
col_order = as.vector(t(outer(c("Plate-based", "Droplet-based"), levels(rank_summary$tier),
                              paste, sep = " | ")))
mat_cons  = reshape2::acast(rank_summary, method ~ column, value.var = "consensus_rank")[, col_order]
mat_mean  = reshape2::acast(rank_summary, method ~ column, value.var = "mean_rank")[, col_order]
mat_worst = reshape2::acast(rank_summary, method ~ column, value.var = "worst_rank")[, col_order]
mat_retr  = reshape2::acast(rank_summary, method ~ column, value.var = "n_retrieved")[, col_order]

# order methods by average consensus rank (not retrieved = last place)
mat_order    = mat_cons
mat_order[is.na(mat_order)] = n_methods
method_order = names(sort(rowMeans(mat_order)))
mat_cons  = mat_cons[method_order, ]
mat_mean  = mat_mean[method_order, ]
mat_worst = mat_worst[method_order, ]
mat_retr  = mat_retr[method_order, ]

# sequential blues; light end stays visibly blue so rank 8 is distinct from "not retrieved" (white)
rank_cols = colorRamp2(c(1, (n_methods + 1) / 2, n_methods), c("#08306B", "#4292C6", "#C6DBEF"))

col_split  = factor(sub(" \\|.*", "", colnames(mat_cons)), levels = c("Plate-based", "Droplet-based"))
col_labels = sub(".*\\| ", "", colnames(mat_cons))

ht_rank = Heatmap(
  mat_cons,
  name = "Consensus rank",
  col = rank_cols,
  na_col = "white",
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split,
  column_labels = col_labels,
  column_names_rot = 45,
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 13),
  rect_gp = gpar(col = "grey85", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (is.na(mat_mean[i, j])) return(invisible(NULL))          # nobody retrieved anything in this column
    if (mat_retr[i, j] == 0) {
      grid.text("–", x, y, gp = gpar(fontsize = 13, col = "grey50"))
    } else {
      txt_col = ifelse(mat_cons[i, j] < (n_methods + 1) / 2, "white", "black")
      grid.text(sprintf("%d\n(%g)", as.integer(mat_cons[i, j]), mat_worst[i, j]), x, y,
                gp = gpar(fontsize = 13, col = txt_col))
    }
  },
  heatmap_legend_param = list(title = "Consensus rank (1 = best)",
                              at = seq_len(n_methods),
                              color_bar = "discrete",
                              direction = "horizontal", nrow = 1,
                              title_position = "leftcenter",
                              labels_gp = gpar(fontsize = 13),
                              title_gp = gpar(fontsize = 13, fontface = "bold"),
                              grid_height = unit(6, "mm"),
                              grid_width  = unit(6, "mm"))
)

png(file.path(path_results_dir, "method_consensus_rank_heatmap.png"),
    width = 14, height = 8, units = "in", res = 300)
draw(ht_rank, heatmap_legend_side = "bottom",
     column_title_gp = gpar(fontsize = 14))
dev.off()
