#' Semi-simulation framework
#'
#' Artificially inflates expression of selected ligand and receptor genes using
#' Negative Binomial sampling, so that both the percentage of cells expressing
#' a gene (PCE) and the magnitude of expression (FC) can be controlled
#' independently.
#'
#' Two design choices worth flagging explicitly:
#'  1. Detection is guaranteed with a "+1" floor for the PCE-selected cells,
#'     and FC-scaled signal is layered on top of that floor. This decouples
#'     "does this cell express the gene at all" (PCE) from "how much" (FC).
#'  2. mu_original must be computed on the count scale (rate * library size),
#'     not on the rate scale alone -- gene_mean from edgeR is a per-library-size
#'     rate (library size already divided out during model fitting), so it
#'     must be multiplied by exp(offset) to get an actual expected count.
#'
#' @param counts     Raw count matrix (genes x cells)
#' @param simulated_interactions_lst  List with $ligand and $receptor gene vectors
#' @param genemetadata  List with $mean (per-CT rates), $disp (dispersion), $edgeR_offsets
#' @param metadata   Cell metadata with cell_ID and Celltype columns
#' @param FC         Fold change multiplier applied to the baseline mean
#' @param PCE        Target percentage of cells (per cell type) to express the gene
#'
#' @return List with:
#'   counts_inflated   - the modified count matrix
#'   PCE_lst           - realized (observed) PCE per gene, for sanity-checking
#'   not_inflated_cells- cells NOT selected for inflation, per cell type
#'   baseline_rate_lst - pre-inflation expression RATE per gene (not counts!),
#'                        used later as the baseline for effective-FC calculations
semi_simulate = function(counts, simulated_interactions_lst, genemetadata, metadata, FC, PCE) {
  
  means_perCT = genemetadata$mean
  offset = genemetadata$edgeR_offsets$edgeR_offset
  
  L_sample = simulated_interactions_lst[["ligand"]]   %>% str_split("_") %>% unlist() %>% unique()
  R_sample = simulated_interactions_lst[["receptor"]] %>% str_split("_") %>% unlist() %>% unique()
  
  counts_inflated    = counts
  PCE_lst             = list()
  not_inflated_cells  = list()
  baseline_rate_lst   = list()   # rate scale -- comparable directly to means_real later
  
  for (tmp_CT in c("Sender", "Receiver")) {
    
    genes_to_sample = if (tmp_CT == "Sender") L_sample else R_sample
    
    # --- Select which cells get inflated (controls PCE) ---
    CT_cells = colnames(counts)[metadata$Celltype == tmp_CT]
    n_cells_to_impute = ceiling(length(CT_cells) * PCE / 100)
    cells_to_impute = sample(CT_cells, n_cells_to_impute)
    not_inflated_cells[[tmp_CT]] = setdiff(CT_cells, cells_to_impute)
    
    cell_offsets = offset[match(cells_to_impute, genemetadata$edgeR_offsets$cell_ID)]
    
    for (gene_sample in genes_to_sample) {
      
      # Wipe this gene's expression across the whole cell type first, so that
      # only the explicitly selected cells end up expressing it (PCE control)
      counts_inflated[gene_sample, CT_cells] = 0
      
      # --- Baseline parameters (both on rate scale, straight from edgeR fit) ---
      gene_rate = means_perCT[means_perCT$gene_names == gene_sample, tmp_CT] %>% as.numeric()
      gene_dispersion = genemetadata$disp %>%
        subset(gene == gene_sample) %>%
        pull(edgeR_dispersion) %>% as.numeric()
      
      # Store the RATE (not count-scale mu) as the pre-inflation baseline.
      # This is what later gets compared against means_real (also rate scale)
      # to compute effective FC
      baseline_rate_lst[[tmp_CT]][[gene_sample]] = gene_rate
      
      # --- Convert rate -> expected count scale for this specific set of cells ---
      # (rate alone can't be sampled from directly -- NB sampling needs real counts)
      mu_original = gene_rate * exp(cell_offsets)
      mu_target   = mu_original * FC
      
      # --- Draw FC-scaled signal, then guarantee detection with a floor of 1 ---
      extra_signal = rnbinom(n_cells_to_impute, mu = mu_target, size = 1 / gene_dispersion)
      counts_inflated[gene_sample, cells_to_impute] = 1 + extra_signal
      
      # Record realized (observed) PCE as a sanity check against the target
      PCE_lst[[tmp_CT]][[gene_sample]] =
        mean(counts_inflated[gene_sample, CT_cells] > 0)
    }
  }
  
  list(
    counts_inflated    = counts_inflated,
    PCE_lst             = PCE_lst,
    not_inflated_cells  = not_inflated_cells,
    baseline_rate_lst   = baseline_rate_lst
  )
}

#' Find the strictest LR-pair rate threshold that still yields at least
#' `target_n` valid pairs
#'
#' Different sequencing technologies vary hugely in sparsity, so a fixed rate
#' threshold doesn't give a comparable number of usable LR pairs across
#' datasets. This scans a range of thresholds and picks the highest one
#' (i.e. best-quality, most-expressed genes) that still clears target_n pairs
#' -- equivalent to finding the n_pairs value closest to target_n from above
#'
#' A pair is "valid" if every ligand gene clears the threshold in Sender
#' cells, and every receptor subunit clears it in Receiver cells (a receptor
#' complex is only functionally detectable if all its subunits are present).
#'
#' @param LR_database  Data frame with "ligand" and "receptor" columns (genes joined by "_")
#' @param means_perCT  Data frame of per-cell-type expression rates, with a "gene_names" column
#' @param offset       Vector of per-cell log-library-size offsets (from edgeR::getOffset)
#' @param metadata     Data frame with a "Celltype" column, aligned to `offset`
#' @param target_n     Minimum number of LR pairs required
#' @param thresholds_to_try  Candidate thresholds to scan, ordered ascending
#'
#' @return A list with: threshold (chosen cutoff), n_pairs (pairs surviving at that
#'   cutoff), valid_rows (logical vector into LR_database), and all_results (the full
#'   threshold -> n_pairs table, useful for diagnostics/plotting)
find_threshold_for_target_pairs = function(counts, LR_database, means_perCT, offset, metadata,
                                           target_n, thresholds_to_try = seq(0.005, 0.5, by = 0.005)) {
  
  # Population-level average library size per cell type.
  mean_libsize_sender   = mean(exp(offset[metadata$Celltype == "Sender"]))
  mean_libsize_receiver = mean(exp(offset[metadata$Celltype == "Receiver"]))
  
  # Check for a given threshold: TRUE for LR pairs where every
  # ligand gene (in Sender) and every receptor gene (in Receiver) has an expected
  # count (rate * mean library size) at or above the threshold.
  check_validity = function(t) {
    apply(LR_database, 1, function(row) {
      ligand_genes   = unlist(str_split(row["ligand"], "_"))
      receptor_genes = unlist(str_split(row["receptor"], "_"))
      all_genes = c(ligand_genes, receptor_genes)
      
      # Skip pairs with genes missing from the dataset entirely
      genes_present = all(all_genes %in% rownames(counts) & all_genes %in% means_perCT$gene_names)
      if (!genes_present) return(FALSE)
      
      ligand_rates   = means_perCT[means_perCT$gene_names %in% ligand_genes,   "Sender"]
      receptor_rates = means_perCT[means_perCT$gene_names %in% receptor_genes, "Receiver"]
      
      # Convert rate -> expected count scale before comparing to threshold
      ligand_mu   = ligand_rates   * mean_libsize_sender
      receptor_mu = receptor_rates * mean_libsize_receiver
      
      all(ligand_mu >= t) & all(receptor_mu >= t)
    })
  }
  
  # Scan all candidate thresholds and record how many pairs survive at each
  results = lapply(thresholds_to_try, function(t) {
    data.frame(threshold = t, n_pairs = sum(check_validity(t)))
  }) %>% do.call(rbind, .)
  
  # Only thresholds that still meet the minimum pair count are usable
  results_ok = results[results$n_pairs >= target_n, ]
  
  if (nrow(results_ok) == 0) {
    stop("No threshold found yielding enough pairs even at the lowest tested value.")
  }
  
  # Pick the threshold whose pair count is closest to target_n.
  # Because n_pairs decreases as threshold increases, this is the strictest
  # (highest) threshold that still clears target_n -- i.e. the best-quality
  # gene set
  best_row = results_ok[which.min(abs(results_ok$n_pairs - target_n)), ]
  best_t   = best_row$threshold
  
  list(
    threshold   = best_t,
    n_pairs     = best_row$n_pairs,
    valid_rows  = check_validity(best_t),
    all_results = results
  )
}

#' Generate AveLogCPM and PCE diagnostic plots to verify simulation signal
#'
#' For a set of benchmark FC/PCE parameter combinations, produces two families
#' of diagnostic plots:
#'   A) AveLogCPM scatter (original vs. inflated), to visually confirm that
#'      simulated ligand/receptor genes show a magnitude shift relative to
#'      the untouched background.
#'   B) A consistency check comparing realized PCE between the lowest and
#'      highest tested FC at a fixed PCE target -- since PCE and FC are meant
#'      to be independent knobs, realized PCE should stay stable across FC
#'      (points should fall near the diagonal).
#'
#' @param counts     Original (pre-simulation) raw count matrix
#' @param master_lst Nested list accumulated across the simulation grid;
#'                    expects master_lst$PCE, $Sender, $Receiver, $counts_aveLogCPM
#'                    keyed by "FC_<fc>_PCE_<pce>" naming
#' @param FC_param   FC values to include in this diagnostic pass (usually min/median/max)
#' @param PCE_param  PCE values to include in this diagnostic pass (usually min/median/max)
#' @param dataset    Dataset name, used only for plot titles
#' @param metadata   Cell metadata with Celltype column
#' @param CT_toPlot  Cell type(s) to include (typically c("Sender", "Receiver"))
#'
#' @return List with avelogcpm_fixedPCE (plots) and PCE_fixedPCE (plots)
compute_diagnostic_plots = function(counts, master_lst, FC_param, PCE_param, dataset, metadata, CT_toPlot) {
  
  plot_avelogcpm_fixed_PCE = list()
  plot_corr_fixed_PCE_cells_expressing = list()
  
  # ---------------------------------------------------------------------
  # 1. Baseline: original AveLogCPM for the target cell types
  # ---------------------------------------------------------------------
  target_cells = metadata %>% filter(Celltype %in% CT_toPlot) %>% pull(cell_ID)
  original_counts_subset = counts[, target_cells]
  original_aveLogCPM = aveLogCPM(original_counts_subset)
  
  # ---------------------------------------------------------------------
  # 2. Loop over each PCE value we want to diagnose
  # ---------------------------------------------------------------------
  for (PCE in PCE_param) {
    
    # master_lst$PCE is keyed by names like "FC_1_PCE_20" -- split on "_" and
    # match the PCE token (4th element) and FC token (2nd element) against
    # the values we're interested in for this diagnostic pass
    param_names = names(master_lst[["PCE"]])
    fc_tokens  = param_names %>% str_split("_") %>% lapply(`[[`, 2)
    pce_tokens = param_names %>% str_split("_") %>% lapply(`[[`, 4)
    
    matches_this_pce = grep(paste0("^", PCE, "$"), pce_tokens)
    matches_requested_fc = which(fc_tokens %in% FC_param)
    keep_idx = intersect(matches_this_pce, matches_requested_fc)
    
    # Identify the lowest- and highest-FC entries among the matches, to use
    # for the PCE consistency scatter plot below
    fc_values_kept = fc_tokens[keep_idx] %>% unlist() %>% as.numeric()
    min_fc_idx = keep_idx[which.min(fc_values_kept)]
    max_fc_idx = keep_idx[which.max(fc_values_kept)]
    
    # --- A. AveLogCPM scatter: original vs. inflated, one plot per FC ---
    for (param_key in param_names[keep_idx]) {
      
      # Ligand + receptor genes get highlighted on the plot, if both cell
      # types are being visualized together
      genes_to_highlight = if (length(CT_toPlot) == 2) {
        c(master_lst$Sender[[param_key]]$L, master_lst$Receiver[[param_key]]$R) %>% unique()
      } else character(0)
      
      inflated_aveLogCPM = master_lst[["counts_aveLogCPM"]][[param_key]]
      current_FC = str_split(param_key, "_")[[1]][2]
      
      df_plot = data.frame(
        original_counts = original_aveLogCPM,
        avelogcpm        = inflated_aveLogCPM,
        is_LR            = names(inflated_aveLogCPM) %in% genes_to_highlight
      )
      
      plot_avelogcpm_fixed_PCE[[param_key]] =
        ggplot(df_plot, aes(x = original_counts, y = avelogcpm, color = is_LR)) +
        geom_point(size = 0.5) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        ggtitle(paste0("PCE = ", PCE, " | Dataset = ", dataset)) +
        xlab("aveLogCPM (Original)") +
        ylab(paste0("aveLogCPM (FC =", current_FC, ")")) +
        theme_light() +
        theme(
          plot.title  = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text  = element_text(size = 12)
        )
    }
    
    # --- B. PCE consistency check: realized PCE at min FC vs. max FC ---
    # If PCE and FC are truly independent, these should sit near the
    # diagonal -- FC should not accidentally be dragging PCE along with it
    df_pce = data.frame(
      FC_Min = master_lst[["PCE"]][[min_fc_idx]],
      FC_Max = master_lst[["PCE"]][[max_fc_idx]]
    )
    
    min_fc_label = str_split(param_names[min_fc_idx], "_")[[1]][2]
    max_fc_label = str_split(param_names[max_fc_idx], "_")[[1]][2]
    
    plot_corr_fixed_PCE_cells_expressing[[paste0("PCE_", PCE)]] =
      ggplot(df_pce, aes(x = FC_Min, y = FC_Max)) +
      geom_point(size = 0.75) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
      ggtitle(paste("PCE Consistency Check | Target PCE =", PCE, "| Dataset:", dataset)) +
      xlab(paste("Realized PCE (FC =", min_fc_label, ")")) +
      ylab(paste("Realized PCE (FC =", max_fc_label, ")")) +
      theme_light() +
      theme(plot.title = element_text(size = 10))
  }
  
  list(
    avelogcpm_fixedPCE = plot_avelogcpm_fixed_PCE,
    PCE_fixedPCE        = plot_corr_fixed_PCE_cells_expressing
  )
}

#' Compare how much of a cell's total sequencing signal comes from the
#' simulated ligand/receptor genes, across the FC/PCE parameter grid
#'
#' For each FC/PCE combination already simulated (reads previously-saved
#' inflated count files from disk), computes what fraction of total counts
#' in the dataset originates from the genes that were artificially inflated.
#'
#' @param dataset   Dataset name (used to locate output files)
#' @param FC_values  FC values to iterate over
#' @param PCE_values PCE values to iterate over
#' @param benchmark_simulated_genes  Character vector of all genes that were
#'   inflated during simulation (ligands + receptor subunits, deduplicated)
#'
#' @return List with $ratio_totalCounts, nested as [[PCE]][[FC]] = percentage
compare_expressionProfile = function(dataset, FC_values, PCE_values, benchmark_simulated_genes) {
  
  sim_dir = file.path("output", paste0(dataset, "_semiSimulation_NB"))
  inflated_counts_files = list.files(sim_dir, pattern = "sc_inflated_counts")
  metadata_files        = list.files(sim_dir, pattern = "sc_metadata")
  
  out_lst = list()
  params_grid = expand.grid(FC = FC_values, PCE = PCE_values)
  
  for (i in 1:nrow(params_grid)) {
    FC  = params_grid[i, "FC"]
    PCE = params_grid[i, "PCE"]
    
    # File naming convention: "..._FC_.*PCE.tsv" -- match the right combination
    file_pattern = paste0("_", FC, "_.*", PCE, ".tsv$")
    inflated_counts_file = inflated_counts_files[grepl(file_pattern, inflated_counts_files)]
    metadata_file         = metadata_files[grepl(file_pattern, metadata_files)]
    
    counts_benchmark = read.table(file.path(sim_dir, inflated_counts_file))
    metadata_benchmark = read.table(file.path(sim_dir, metadata_file))
    
    # What fraction of all sequenced counts came from the simulated LR genes?
    is_LR_gene = rownames(counts_benchmark) %in% benchmark_simulated_genes
    total_counts_from_LR = sum(counts_benchmark[is_LR_gene, ])
    total_library_size   = sum(counts_benchmark)
    
    magnitude_pct = if (total_library_size > 0) {
      (total_counts_from_LR / total_library_size) * 100
    } else {
      0
    }
    
    out_lst[["ratio_totalCounts"]][[paste0("PCE", PCE)]][[paste0("FC", FC)]] = magnitude_pct
  }
  
  gc()
  
  out_lst
}


###################################################
### Helper functions for Visualization.R script ###
###################################################

#' Load the simulated-interactions ground truth, inflated counts, and realized
#' PCE for one FC/PCE combination
load_param_combination = function(file_path, inflated_counts_files, PCE_files,
                                  simulated_interactions, FC_val, PCE_val) {
  
  param_naming = paste0("FC_", FC_val, "_PCE_", PCE_val)
  file_pattern = paste0("_", FC_val, "_.*", PCE_val)
  
  counts_file = inflated_counts_files[grepl(paste0(file_pattern, ".tsv$"), inflated_counts_files)]
  inflated_counts = read.table(file.path(file_path, counts_file), header = TRUE, row.names = 1)
  
  pce_file = PCE_files[grepl(paste0(file_pattern, ".RDS$"), PCE_files)]
  realized_pce_pct = readRDS(file.path(file_path, pce_file)) %>% unlist() * 100
  
  list(
    param_naming     = param_naming,
    inflated_counts  = inflated_counts,
    ligand_receptor  = simulated_interactions[["ligand_receptor"]],
    realized_pce_pct = realized_pce_pct
  )
}

#' Split "Ligand_Receptor1_Receptor2_Receptor3" interaction strings into
#' component genes and return the unique ligand/receptor gene sets
extract_inflated_genes = function(ligand_receptor_strings) {
  
  split_table = ligand_receptor_strings %>%
    str_split_fixed("_", n = 4) %>%
    as.data.frame()
  colnames(split_table) = c("ligand", "receptor1", "receptor2", "receptor3")
  
  list(
    ligands   = split_table$ligand %>% unique() %>% setdiff(""),
    receptors = c(split_table$receptor1, split_table$receptor2, split_table$receptor3) %>%
      unique() %>% setdiff("")
  )
}

#' Sample n_genes from a cell type's inflated gene set and build density plots for each
generate_gene_density_plots = function(genes_for_CT, CT, inflated_counts, original_counts,
                                       metadata, FC_val, PCE_val, n_genes = 9) {
  
  stopifnot(metadata$cell_ID == colnames(inflated_counts))
  
  set.seed(1)  # reproducible gene sampling across runs
  sampled_genes = genes_for_CT %>% unlist() %>% sample(n_genes)
  
  is_sender = (CT == "Sender")
  
  setNames(
    lapply(sampled_genes, build_gene_density_plot,
           inflated_counts = inflated_counts, original_counts = original_counts,
           metadata = metadata, is_sender = is_sender, FC_val = FC_val, PCE_val = PCE_val),
    sampled_genes
  )
}

#' Build a single count-density diagnostic plot (inflated vs. original) for one gene
build_gene_density_plot = function(gene, inflated_counts, original_counts, metadata,
                                   is_sender, FC_val, PCE_val) {
  
  cell_mask = if (is_sender) (metadata$Celltype == "Sender") else (metadata$Celltype == "Receiver")
  
  tmp_df = data.frame(
    inflated_counts = as.numeric(inflated_counts[gene, cell_mask]),
    original_counts = as.numeric(original_counts[gene, cell_mask])
  ) %>% reshape2::melt()
  
  group_means = plyr::ddply(tmp_df, "variable", plyr::summarise, grp.mean = mean(value))
  
  p = ggplot(tmp_df, aes(x = value, color = variable)) +
    geom_density() +
    geom_vline(data = group_means, aes(xintercept = grp.mean, color = variable), linetype = "dashed") +
    ggtitle(paste0("Count density | FC:", FC_val, " PCE:", PCE_val, " Gene:", gene)) +
    theme_light() +
    theme(
      plot.title   = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.text  = element_text(size = 12),
      axis.text    = element_text(size = 10)
    )
  
  p$plot_env = rlang::base_env()   # remvoe env to save memory
  p
}

#' Build the dataset-level mean-variance diagnostic plot, highlighting genes
#' sampled as baseline vs. selected as simulated ligands/receptors
build_meanvar_plot = function(genemetadata, simulated_ligands, simulated_receptors, dataset) {
  
  tmp_df_L = subset(genemetadata$mean, gene_names %in% simulated_ligands)
  tmp_df_R = subset(genemetadata$mean, gene_names %in% simulated_receptors)
  baseline_df = genemetadata$mean[genemetadata$mean$randomly_sampled, ]
  
  ggplot(genemetadata$mean, aes(x = means, y = vars)) +
    geom_hex() +
    # Red: randomly sampled baseline (10 genes per expression bin)
    geom_point(data = baseline_df, aes(x = means, y = vars), colour = "red", size = 2) +
    geom_smooth(data = baseline_df, aes(x = means, y = vars), method = "lm", color = "red") +
    # Green: simulated ligands
    geom_point(data = tmp_df_L, aes(x = means, y = vars), colour = "green", size = 2) +
    geom_smooth(data = tmp_df_L, aes(x = means, y = vars), method = "lm", color = "green") +
    # Yellow: simulated receptors
    geom_point(data = tmp_df_R, aes(x = means, y = vars), colour = "yellow", size = 2) +
    geom_smooth(data = tmp_df_R, aes(x = means, y = vars), method = "lm", color = "yellow") +
    labs(
      title = paste0("Mean-Variance Profile | Dataset: ", dataset),
      subtitle = "Red: Random | Green: Ligands | Yellow: Receptors",
      x = "Mean (log10)", y = "Variance (log10)"
    ) +
    theme_light()
}

#' Combine per-dataset, per-method metric results into one long dataframe,
#' attaching the realized median FC (from the diagnostic FC step) to each row
compute_recallprecision_stats = function(metric_results, diagnostic_df_realFC) {
  stats = lapply(names(metric_results), function(dataset) {
    tmp_diagnostic_df_realFC = diagnostic_df_realFC[[dataset]]
    
    lapply(metric_results[[dataset]], function(method) {
      method %>% do.call(rbind, .) %>%
        left_join(tmp_diagnostic_df_realFC %>% select(FC, PCE, FC_real_median), by = c("FC", "PCE")) %>%
        mutate(dataset = dataset)
    }) %>% do.call(rbind.data.frame, .)
  })
  names(stats) = names(metric_results)
  stats
}

#' Extract realized effective-FC values for one FC/PCE combination from the
#' saved realFC_aftersimulation RDS file
#'
#' @param file_path Directory containing this dataset's realFC_aftersimulation_*.RDS files
#' @param FC_val    Theoretical FC value for this combination
#' @param PCE_val   Theoretical PCE value for this combination
#'
#' @return Data frame with columns: FC, PCE, gene, effective_FC
extract_effective_fc = function(file_path, FC_val, PCE_val) {
  
  realFC_file = list.files(file_path, pattern = "realFC_aftersimulation") %>%
    magrittr::extract(grepl(paste0("_", FC_val, "_.*", PCE_val, ".RDS$"), .))
  
  realized_fc = readRDS(file.path(file_path, realFC_file)) %>% unlist() %>% subset(!is.infinite(.))
  
  data.frame(
    FC           = FC_val,
    PCE          = PCE_val,
    gene         = names(realized_fc),
    effective_FC = as.numeric(realized_fc)
  )
}

#' Build the collapsed effective-FC-vs-theoretical-FC boxplot (faceted by PCE)
#' across the full FC/PCE grid, plus the per-combination summary stats needed
#' downstream by compute_recallprecision_stats()
#'
#' @param df       Combined output of extract_effective_fc() across the whole grid
#'                 (rbind of one call per FC/PCE combination)
#' @param dataset  Dataset name, used only for the plot title
#'
#' @return List with: plot (single ggplot, faceted by PCE) and summary
#'   (data.frame with FC, PCE, theoreticalFC, PCE_real, FC_real_median --
#'   same shape previously produced per-combination and rbind'ed manually)
plot_FCafter_semisimulation = function(df, dataset) {
  
  fc_levels = sort(unique(df$FC))
  df$FC_factor = factor(df$FC, levels = fc_levels)
  df$PCE_label = factor(paste0("PCE = ", df$PCE), levels = paste0("PCE = ", sort(unique(df$PCE))))
  
  # One reference point per FC level, marking the theoretical target
  ref_df = data.frame(FC_factor = factor(fc_levels, levels = fc_levels), FC = fc_levels)
  
  plot = ggplot(df, aes(x = FC_factor, y = effective_FC)) +
    geom_boxplot(outlier.shape = NA, fill = "grey95", width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.7) +
    geom_point(data = ref_df, aes(x = FC_factor, y = FC),
               color = "red", shape = 95, size = 8, inherit.aes = FALSE) +
    geom_line(data = ref_df, aes(x = FC_factor, y = FC, group = 1),
              color = "red", linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
    facet_wrap(~ PCE_label) +
    coord_trans(y = "log10") +   # effective FC can span well below 1 to very large outliers
    theme_light() +
    labs(
      title = paste0("Effective vs. theoretical fold change | Dataset: ", dataset),
      x = "Theoretical FC", y = "Effective FC (log scale)"
    ) +
    theme(
      plot.title  = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      strip.text  = element_text(size = 12)
    )
  
  summary = df %>%
    group_by(FC, PCE) %>%
    dplyr::summarise(theoreticalFC = first(FC), PCE_real = first(PCE),
                     FC_real_median = median(effective_FC), .groups = "drop") %>%
    as.data.frame()
  
  list(plot = plot, summary = summary)
}

#' Read one method's significant-interactions output file and return the
#' ligand_receptor IDs flagged as significant
read_significant_ligand_receptors = function(filepath) {
  read.table(filepath, header = TRUE) %>%
    mutate(significant = as.logical(significant)) %>%
    filter(significant) %>%
    pull(ligand_receptor) %>%
    unname()
}

#' Build a method x interaction binary membership matrix: for each method,
#' collapse across all its files (e.g. FC/PCE combinations) into "did this
#' method ever call this interaction significant" (1/0), for every ID in
#' `universe`.
build_method_membership_matrix = function(interactions_by_method_and_file, universe) {
  lapply(interactions_by_method_and_file, function(files_for_method) {
    per_file_binary = as.data.frame(lapply(files_for_method, function(x) as.numeric(universe %in% x)))
    data.frame(method = ifelse(rowSums(per_file_binary) > 0, 1, 0), row.names = universe)
  }) %>% do.call(cbind, .)
}

#' Build FC/PCE column annotation metadata (factor levels + color palettes)
#' from heatmap column names like "FC:1\nPCE:20" or "FC_1-PCE_20"
build_fc_pce_col_meta = function(colnames_vec, sep, clean_fn) {
  col_meta = data.frame(colnames_vec) %>%
    separate(1, into = c("FC", "PCE"), sep = sep) %>%
    dplyr::mutate(across(everything(), clean_fn))
  
  col_meta$FC  = factor(col_meta$FC,  levels = unique(col_meta$FC))
  col_meta$PCE = factor(col_meta$PCE, levels = unique(col_meta$PCE))
  
  list(
    FC  = col_meta$FC,
    PCE = col_meta$PCE,
    fc_cols  = setNames(brewer.pal(length(unique(col_meta$FC)),  "Set1"), unique(col_meta$FC)),
    pce_cols = setNames(brewer.pal(length(unique(col_meta$PCE)), "Set2"), unique(col_meta$PCE))
  )
}

#' F1-score heatmap: methods x FC/PCE combinations
build_f1_heatmap = function(stats_df) {
  mat = stats_df %>%
    mutate(Param_Comb = paste0("FC:", FC, "\nPCE:", PCE)) %>%
    select(method, Param_Comb, f1score) %>%
    pivot_wider(names_from = Param_Comb, values_from = f1score) %>%
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  col_meta = build_fc_pce_col_meta(colnames(mat), sep = "\n", clean_fn = ~gsub("FC:", "", .))
  top_ann  = HeatmapAnnotation(FC = col_meta$FC, col = list(FC = col_meta$fc_cols, PCE = col_meta$pce_cols))
  col_f1   = colorRamp2(c(0, 0.20, 0.60, 1), c("grey95", "#7fcdbb", "#2c7fb8", "#081d58"))
  
  heatmap = Heatmap(mat,
                    name = "F1 Score", col = col_f1, top_annotation = top_ann,
                    column_split = col_meta$PCE, cluster_rows = FALSE, cluster_columns = FALSE,
                    width = ncol(mat) * unit(6, "mm"), height = nrow(mat) * unit(6, "mm"),
                    row_names_side = "left", column_names_gp = gpar(fontsize = 8),
                    show_column_names = FALSE, rect_gp = gpar(col = "white", lwd = 1))
  
  list(heatmap = heatmap, mat = mat)
}

# Detection completeness versus ranking quality across CCC methods
build_detection_vs_rank_plot = function(stats_df) {
  ggplot(stats_df, aes(x = detection_rate, y = normalizedRank_conditional, color = method)) +
    geom_point(size = 2, alpha = 0.8) +
    facet_grid(~PCE) +
    scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
    labs(
      title    = "Detection completeness vs. ranking quality by method",
      subtitle = "Lower-right = detects most simulated interactions and ranks them highly",
      x = "Detection rate (% of simulated interactions reported)",
      y = "Mean normalized rank (conditional on detection, 0 = best)"
    ) +
    theme_light()
}

#' Sum, across every pair of methods, how many interactions both methods
#' flagged (both columns == 1) in a binary membership dataframe
get_pairwise_intersects = function(df) {
  method_names = colnames(df)
  pairs = combn(method_names, 2)
  
  results = apply(pairs, 2, function(p) {
    sum(df[[p[1]]] == 1 & df[[p[2]]] == 1)
  })
  names(results) = apply(pairs, 2, paste, collapse = " & ")
  results
}

#' Heatmap of % of simulated interactions retrieved (pairwise-intersection
#' summary), across the FC/PCE grid. `ref_dims` reuses the F1 heatmap's
#' matrix dimensions for consistent cell sizing (matches original coupling).
build_retrieval_pct_heatmap = function(lst_binarydf, total_simulated, ref_dims) {
  summary_mat = lapply(lst_binarydf, get_pairwise_intersects) %>% do.call(cbind, .)
  summary_mat_pct = (summary_mat / total_simulated) * 100
  
  colnames(summary_mat_pct) = colnames(summary_mat_pct) %>% gsub("FC_", "", .)
  col_meta = build_fc_pce_col_meta(colnames(summary_mat_pct), sep = "-", clean_fn = ~gsub("_", ":", .))
  top_ann  = HeatmapAnnotation(FC = col_meta$FC, col = list(FC = col_meta$fc_cols, PCE = col_meta$pce_cols))
  col_f1   = colorRamp2(c(0, 20, 60, 100), c("grey95", "#7fcdbb", "#2c7fb8", "#081d58"))
  
  Heatmap(summary_mat_pct,
          name = "% significant int. retrieved", col = col_f1,
          cluster_rows = FALSE, top_annotation = top_ann, column_split = col_meta$PCE,
          cluster_columns = FALSE, row_names_side = "left",
          width = ref_dims$ncol * unit(6, "mm"), height = ref_dims$nrow * unit(20, "mm"),
          column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 12),
          show_column_names = FALSE, rect_gp = gpar(col = "white", lwd = 1))
}

#' Dot-plot alternative to the F1 heatmap: one row per method, faceted by PCE
build_f1_dotplot = function(stats_df) {
  method_order = stats_df %>%
    group_by(method) %>%
    dplyr::summarise(m = mean(f1score)) %>%
    arrange(m) %>%
    pull(method)
  
  df = stats_df
  df$method = factor(df$method, levels = method_order)
  df$PCE_label = factor(df$PCE, levels = c(4, 7, 10, 20, 30, 40))
  
  ggplot(df, aes(x = f1score, y = method)) +
    geom_segment(aes(x = 0, xend = 1, y = method, yend = method), color = "gray90", size = 0.5) +
    geom_point(aes(color = factor(FC)), size = 3, alpha = 0.8) +
    facet_wrap(~ PCE_label, ncol = 2) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(
      panel.spacing = unit(1, "lines"),
      strip.background = element_rect(fill = "gray95", color = NA),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10))
    ) +
    labs(title = "Comparing F1 scores across sparsity (PCE) and signal (FC)",
         x = "F1 Score", y = "", color = "Fold Change")
}

#' Generic UpSet plot builder (shared styling used for all three UpSet plots
#' in this script: per-FC/PCE grid, all-significant, and intersect-simulated)
build_upset_plot = function(binary_df, title) {
  m = make_comb_mat(binary_df)
  UpSet(m,
        top_annotation = upset_top_annotation(m, height = unit(10, "cm"), add_numbers = TRUE,
                                              numbers_gp = gpar(fontsize = 12), gp = gpar(fill = "steelblue")),
        pt_size = unit(5, "mm"),
        lwd = 1.5,
        right_annotation = upset_right_annotation(m, width = unit(4, "cm"), gp = gpar(fill = "darkred")),
        row_names_gp = gpar(fontsize = 12),
        row_gap = unit(0, "mm"),
        column_title = title,
        column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
}

#' Build per-FC/PCE-combination UpSet plots (and their underlying binary
#' membership tables) showing which methods retrieved which simulated
#' interactions
build_grid_upset_data = function(significant_interactions_lst, FC, PCE) {
  
  known_methods = c("cellphonedb", "singlecellsignalR", "cellchat", "connectome",
                    "natmi_specificity", "geometric_mean", "scseqcomm", "log2fc")
  
  lst_upset_plots_FC  = list()
  lst_upset_plots_PCE = list()
  lst_binarydf        = list()
  
  for (fc in FC) {
    for (pce in PCE) {
      
      file_key = paste0("significant_interactions_FC_", fc, "_PercCellsExpressing_", pce, ".csv")
      per_method = lapply(significant_interactions_lst$intersect_significant_simulated,
                          function(method) method[[file_key]])
      
      grid_key = paste0("FC_", fc, "-PCE_", pce)
      
      if (length(unlist(per_method)) == 0) {
        # No method retrieved a single simulated interaction at this FC/PCE
        binary_df = as.data.frame(setNames(rep(list(NA), length(known_methods)), known_methods))
        lst_binarydf[[grid_key]] = binary_df
        lst_upset_plots_FC[[paste0("FC_", fc)]][[paste0("PCE_", pce)]]   = 0
        lst_upset_plots_PCE[[paste0("PCE_", pce)]][[paste0("FC_", fc)]] = 0
      } else {
        all_strings = unique(unlist(per_method))
        binary_df = as.data.frame(lapply(per_method, function(x) as.numeric(all_strings %in% x)))
        rownames(binary_df) = all_strings
        lst_binarydf[[grid_key]] = binary_df
        
        plot = build_upset_plot(binary_df, paste0("All significant interactions retrieved for FC ", fc, " PCE ", pce))
        lst_upset_plots_FC[[paste0("FC_", fc)]][[paste0("PCE_", pce)]]   = plot
        lst_upset_plots_PCE[[paste0("PCE_", pce)]][[paste0("FC_", fc)]] = plot
      }
    }
  }
  
  list(lst_upset_plots_FC = lst_upset_plots_FC, lst_upset_plots_PCE = lst_upset_plots_PCE, lst_binarydf = lst_binarydf)
}

#' Build one precision-recall plot per unique value of a facet column (PCE or FC)
build_precision_recall_plots = function(df, facet_var) {
  # unique() (not sort) preserves first-appearance order, matching the
  # original loop's list-key insertion order (which iterated the raw column
  # -- including duplicate rows -- and simply overwrote repeated keys)
  facet_values = unique(df[[facet_var]])
  
  plots = lapply(facet_values, function(val) {
    data = df[df[[facet_var]] == val, ]
    ggplot(data, aes(y = precision, x = recall, color = method)) +
      geom_point(size = 1.5) +
      geom_line() +
      theme_light() +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
      theme(legend.title = element_text(size = 7), legend.text = element_text(size = 7),
            axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
            plot.title = element_text(size = 10)) +
      ggtitle(paste0(facet_var, "_", val))
  })
  
  setNames(plots, paste0(facet_var, "_", facet_values))
}

#' F1-score-vs-parameter and normalized-rank-vs-FC plots
build_f1_and_rank_plots = function(stats_df) {
  p1 = ggplot(stats_df, aes(x = PCE, y = f1score, color = method)) +
    geom_point() + geom_line() + facet_grid(~FC) + ggtitle("Faceted by FC") + theme_light() +
    theme(
      plot.title   = element_text(size = 14),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.text  = element_text(size = 14),
      legend.title  = element_text(size = 14),
      axis.text    = element_text(size = 12),
      strip.text = element_text(size = 14)
    )
  
  p2 = ggplot(stats_df, aes(x = FC, y = f1score, color = method)) +
    geom_point() + geom_line() + facet_grid(~PCE) + ggtitle("Faceted by PCE") + theme_light()+
    theme(
      plot.title   = element_text(size = 14),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.text  = element_text(size = 14),
      legend.title  = element_text(size = 14),
      axis.text    = element_text(size = 12),
      strip.text = element_text(size = 14)
    )
  
  p3 = ggplot(stats_df, aes(x = FC, y = normalizedRank_simulated_among_significant, color = method)) +
    geom_point() + geom_line() + facet_grid(~PCE) + ggtitle("Faceted by PCE") + theme_light()+
    ylab("Mean normalized rank of simulated interactions") +
    theme(
      plot.title   = element_text(size = 14),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.text  = element_text(size = 14),
      legend.title  = element_text(size = 14),
      axis.text    = element_text(size = 12),
      strip.text = element_text(size = 14)
    )
  
  list(f1_by_FC = p1, f1_by_PCE = p2, rank_by_PCE = p3)
}

#' Build the relative-F1-score matrix used for the cross-dataset stacked heatmap
build_relative_f1_matrix = function(statistics_results_lst_recallprecision_plot) {
  x = do.call(rbind.data.frame, statistics_results_lst_recallprecision_plot) %>%
    mutate(PCE = as.numeric(as.character(PCE)), FC = as.numeric(as.character(FC))) %>%
    arrange(PCE, FC) %>%
    mutate(params = paste0("PCE_", PCE, "_FC_", FC)) %>%
    mutate(params = factor(params, levels = unique(params)))
  
  df_relative = x %>%
    group_by(dataset, params) %>%
    dplyr::mutate(relative_score = f1score / sum(f1score, na.rm = TRUE)) %>%
    dplyr::summarise(rel_vec = list(setNames(relative_score, method)), .groups = "drop")
  
  mat_rel = df_relative %>%
    pivot_wider(names_from = params, values_from = rel_vec) %>%
    tibble::column_to_rownames("dataset") %>%
    as.matrix()
  
  list(x = x, mat_rel = mat_rel)
}

#' Stacked-segment heatmap showing each method's relative share of F1 score,
#' per dataset (row) x PCE/FC combination (column)
build_relative_f1_stacked_heatmap = function(mat_rel, method_cols) {
  column_pce_values = str_split_i(colnames(mat_rel), "_", 2) %>% str_c("PCE:", .)
  column_fc_labels  = str_split_i(colnames(mat_rel), "_", 4) %>% str_c("FC:", .)
  column_pce_values = factor(column_pce_values, levels = unique(column_pce_values))
  
  Heatmap(mat_rel,
          show_heatmap_legend = FALSE,
          col = method_cols,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = factor(rownames(mat_rel), levels = rownames(mat_rel)),
          row_title_side = "right",
          row_title_rot = 0,
          row_gap = unit(2, "mm"),
          show_row_names = FALSE,
          column_split = column_pce_values,
          column_title_side = "top",
          column_gap = unit(5, "mm"),
          column_labels = column_fc_labels,
          column_names_side = "bottom",
          column_names_rot = 90,
          column_names_gp = gpar(fontsize = 12),
          border = TRUE,
          rect_gp = gpar(type = "none", col = "black", lwd = 1),
          cell_fun = function(j, i, x, y, width, height, fill) {
            scores = mat_rel[[i, j]]
            if (is.null(scores) || any(is.na(scores))) {
              grid.rect(x, y, width * 0.9, height, gp = gpar(fill = "grey90", col = NA))
              return()
            }
            current_y = y - height / 2
            for (k in seq_along(scores)) {
              seg_h = height * scores[k]
              grid.rect(x = x, y = current_y + seg_h / 2, width = width * 0.9, height = seg_h,
                        gp = gpar(fill = method_cols[names(scores)[k]], col = NA))
              current_y = current_y + seg_h
            }
          })
}
