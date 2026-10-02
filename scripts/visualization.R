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
  #lst_precision_recall_byPCE = build_precision_recall_plots(statistics_results_lst_recallprecision_plot[[dataset]], "PCE")
  #lst_precision_recall_byFC  = build_precision_recall_plots(statistics_results_lst_recallprecision_plot[[dataset]], "FC")
  
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
  
  # --- Save all plots for this dataset ---
  pdf(file.path(path_results_dir, paste0(dataset, "_results_plots.pdf")), width = 15, height = 7)
  #wrap_plots(lst_precision_recall_byPCE, guides = "collect") %>% print
  #wrap_plots(lst_precision_recall_byFC, guides = "collect")  %>% print
  summary_plots$f1_by_FC   %>% print
  summary_plots$f1_by_PCE  %>% print
  summary_plots$rank_by_PCE %>% print
  p1 %>% print
  heatmap1 %>% print
  #heatmap1_alternative %>% print
  heatmap2 %>% print
  
  lst_upset_plots_FC$FC_1$PCE_20  %>% print
  lst_upset_plots_FC$FC_7$PCE_20 %>% print
  
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

pdf(file.path(path_results_dir, "results_acrossdatasets_plots.pdf"), width = 13, height = 10)
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
# 6. Create recommendation hierarchical tree
# =============================================================================

suppressMessages({
  library(dplyr)
  library(stringr)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
})

# Input data (per-dataset FC-threshold -> expression-signal-fraction results)
df_smartseq2 = data.frame(
  dataset    = "SMARTseq2",
  PCE        = c(7,7,20,20,40,40),
  Expression = c("< 0.349%","> 0.349%","< 0.370%","> 0.370%","< 0.421%","> 0.421%") ,
  Method     = c("-","log2fc",
                 "geometric_mean|natmi_specificity|cellphonedb","log2fc|geometric_mean",
                 "cellchat","cellchat")
)
# (fix: keep Expression labels simple and explicit instead of the regex trick above)
df_smartseq2$Expression = c("< 0.349%","> 0.349%","< 0.370%","> 0.370%","< 0.421%","> 0.421%")

df_10x = data.frame(
  dataset    = "10x",
  PCE        = c(7,7,20,20,40,40),
  Expression = c("< 0.468%","> 0.468%","< 0.576%","> 0.576%","< 0.755%","> 0.755%"),
  Method     = c("geometric_mean","geometric_mean",
                 "scseqcomm|singlecellsignalR","scseqcomm|singlecellsignalR",
                 "scseqcomm|singlecellsignalR|cellchat","scseqcomm|cellchat|singlecellsignalR")
)

df_vasaseq = data.frame(
  dataset    = "VASAseq",
  PCE        = c(7,7,20,20,40,40),
  Expression = c("< 1.42%","> 1.42%","< 1.73%","> 1.73%","< 2.23%","> 2.23%"),
  Method     = c("-","log2fc",
                 "log2fc","log2fc",
                 "singlecellsignalR","log2fc|singlecellsignalR|scseqcomm")
)

df_BD = data.frame(
  dataset    = "BDRhapsody",
  PCE        = c(7,7,20,20,40,40),
  Expression = c("< 1.55%","> 1.55%","< 1.79%","> 1.79%","< 2.22%","> 2.22%"),
  Method     = c("geometric_mean","geometric_mean|singlecellsignalR",
                 "singlecellsignalR|geometric_mean|cellphonedb","singlecellsignalR",
                 "singlecellsignalR|cellchat","singlecellsignalR|cellchat|scseqcomm|geometric_mean")
)

df_10x_immune_R1 = data.frame(
  dataset    = "10x_immune_R1",
  PCE        = c(7,7,20,20,40,40),
  Expression = c("< 1.78%","> 1.78%","< 2.04%","> 2.04%","< 2.42%","> 2.42%"),
  Method     = c("singlecellsignalR|geometric_mean","singlecellsignalR",
                 "scseqcomm|singlecellsignalR","scseqcomm|singlecellsignalR",
                 "scseqcomm|cellchat","scseqcomm|cellchat")
)

all_datasets = list(
  SMARTseq2      = df_smartseq2,
  `10x`          = df_10x,
  VASAseq        = df_vasaseq,
  BDRhapsody     = df_BD,
  `10x_immune_R1` = df_10x_immune_R1
)

# Helper functions

#' Check whether two "method1|method2|..." strings represent the same set of
#' methods, ignoring order -- used to collapse a PCE branch into a single
#' "Any Expr" leaf when the winning method(s) don't actually change across
#' the expression threshold
same_method_set = function(m1, m2) {
  setequal(str_split(m1, "\\|")[[1]], str_split(m2, "\\|")[[1]])
}

#' Convert a "method1|method2|..." string into a Graphviz-friendly multi-line
#' label
method_label = function(m) {
  str_replace_all(m, "\\|", "\\\\n")
}

#' Sanitize a dataset name into a valid Graphviz graph/node ID
safe_id = function(x) {
  gsub("[^A-Za-z0-9]", "_", x)
}

# Build the DOT source for one dataset's tree
build_tree_dot = function(dataset_name, df) {
  
  graph_id = paste0("selection_tree_", safe_id(dataset_name))
  pces = sort(unique(df$PCE))
  
  lines = c(
    sprintf("digraph %s {", graph_id),
    # tighter overall spacing between ranks/nodes, smaller margins
    '  graph [layout = dot, ranksep = 0.25, nodesep = 0.15, bgcolor = "white", margin = 0]',
    '  node [shape = box, fontname = Helvetica, style = filled, fillcolor = "white", ',
    '        margin = "0.04,0.02", fixedsize = false]',
    '  edge [color = "#566573", arrowhead = vee, fontname = Helvetica, fontsize = 8, penwidth = 0.8]',
    "",
    sprintf('  Root [label = "%s", fillcolor = "#AED6F1", style = "filled,bold", fontsize = 10, height = 0.25, width = 0.9]', dataset_name),
    ""
  )
  
  # --- PCE nodes (smaller, tighter) ---
  lines = c(lines, '  node [fillcolor = "#FCF3CF", height = 0.22, width = 0.45, fontsize = 9]')
  pce_ids = setNames(paste0("PCE_", pces), pces)
  for (pce in pces) {
    lines = c(lines, sprintf('  %s [label = "%s%%"]', pce_ids[[as.character(pce)]], pce))
  }
  lines = c(lines, "", '  node [fillcolor = "white", fontsize = 7, style = filled, height = 0.2]')
  lines = c(lines, sprintf("  Root -> {%s}", paste(pce_ids, collapse = " ")), "")
  
  # --- Method leaves per PCE tier ---
  leaf_counter = 0
  for (pce in pces) {
    rows = df %>% filter(PCE == pce)
    row_low  = rows[1, ]
    row_high = rows[2, ]
    pid = pce_ids[[as.character(pce)]]
    
    leaf_counter = leaf_counter + 1
    leaf_low_id = paste0("M_", leaf_counter)
    lines = c(lines,
              sprintf('  %s [label = "%s"]', leaf_low_id, method_label(row_low$Method)),
              sprintf('  %s -> %s [label = "%s"]', pid, leaf_low_id, row_low$Expression)
    )
    
    leaf_counter = leaf_counter + 1
    leaf_high_id = paste0("M_", leaf_counter)
    lines = c(lines,
              sprintf('  %s [label = "%s"]', leaf_high_id, method_label(row_high$Method)),
              sprintf('  %s -> %s [label = "%s"]', pid, leaf_high_id, row_high$Expression)
    )
  }
  
  lines = c(lines, "}")
  paste(lines, collapse = "\n")
}

# Render each dataset's tree to its own standalone SVG
output_dir = "output/results/decision_trees"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (dataset_name in names(all_datasets)) {
  df = all_datasets[[dataset_name]]
  
  dot_src = build_tree_dot(dataset_name, df)
  graph = grViz(dot_src)
  
  svg_text = export_svg(graph)
  svg_path = file.path(output_dir, paste0(dataset_name, ".svg"))
  writeLines(svg_text, svg_path)
  
  # Optional: also export a quick PDF/PNG for preview
  rsvg_pdf(charToRaw(svg_text), file = file.path(output_dir, paste0(dataset_name, ".pdf")))
  
  message("Built: ", svg_path)
}

