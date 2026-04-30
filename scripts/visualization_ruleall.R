suppressMessages({
  # Load package
  library(dplyr)
  library(optparse)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  library(purrr) 
  library(magrittr)
  library(plyr)
  library(tidyr)
  library(ComplexUpset)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(patchwork)
  library(reshape2)
  library(edgeR)
  source("scripts/helper_functions.R")
}) 

# Get list with command line arguments by name
option_list = list(
  make_option(c("--config.yaml"), type="character", default=NULL, help="config file", metavar="character"),
  make_option(c("--path_output_dir"), type="character", default=NULL, help="path to output directory", metavar="character"),
  make_option(c("--path_results_dir"), type="character", default=NULL, help="path to where to save the results", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# An useful variability if the argument is missing
if (is.null(opt$config.yaml) | is.null(opt$path_output_dir) | is.null(opt$path_results_dir)){
  print_help(opt_parser)
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Call the argument
path_config.yaml = opt$config.yaml
path_output_dir = opt$path_output_dir
path_results_dir = opt$path_results_dir
dir.create(path_results_dir)

########################################
##### Loading and processing files #####
########################################
metric_results = readRDS("output/final_scores.RDS")

config = yaml::read_yaml(path_config.yaml)

# Setting parameters
methods = config$methods %>% unlist
datasets = config$datasets %>% unlist
FC = config$semiSimulation$FC %>% unlist %>% as.double
PCE = config$semiSimulation$PCE %>% unlist %>% as.integer

############################
##### Diagnostic plots #####
############################

diagnostic_plots_perCT= list()
diagnostic_plots_MeanVar = list()
diagnostic_plots_realFC = list()
diagnostic_df_realFC = list()
diagnostic_gene_densityPlots = list()
simulated_interactions_lst = list() # selected LR genes to inflate counts are all the same across parameters. Only change with dataset
for(dataset in datasets)
{
  # load files
  file_path = file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB"))
  inflated_counts_files = file_path %>% list.files(., pattern = "sc_inflated_counts")
  original_counts = read.table(file.path("data/",dataset, "raw_counts.tsv"))  # load original counts
  rownames(original_counts) = toupper(rownames(original_counts))
  colnames(original_counts) = gsub("[.]", "_" , colnames(original_counts))
  
  PCE_files = file_path %>% list.files(., pattern = "PCE_perGene")
  simulated_interactions_files = file_path %>% list.files(., pattern = "simulated_interactions")
  
  # list to save all data
  master_lst = list()
  
  # load metadata files
  # as the metadata files are the same -> load 1st one
  metadata_path = file.path("data/processed",dataset)
  metadata_file = list.files(metadata_path, pattern = "metadata_")
  metadata = read.table(file.path(metadata_path,metadata_file))
  
  metadata$cell_ID = gsub("[.-]","_" , metadata$cell_ID)
  
  # generate the grid of parameters used for naming the list to generate outputs
  params_grid = expand.grid(vector1 = FC, vector2 = PCE)
  
  # iterate over the grid of parameters
  for(i in 1:nrow(params_grid))
  {
    x = params_grid[i,] %>% as.numeric ; names(x) = c("FC","PCE")
    
    naming = paste0("FC_",x["FC"],"_PCE_",x["PCE"])
    
    # load correct count file depending on the params_grid
    inflated_counts_file = inflated_counts_files %>% magrittr::extract(grepl(paste0("_" , x["FC"] , "_", ".*",x["PCE"] , ".tsv$"), inflated_counts_files))
    inflated_counts = read.table(file.path(file_path,inflated_counts_file))
    
    # load correct simulated interactions file depending on the params_grid
    simulated_interactions_file = simulated_interactions_files %>% magrittr::extract(grepl(paste0("_" , x["FC"] , "_",".*",x["PCE"] , ".RDS$"), simulated_interactions_files))
    tmp_vec = readRDS(file.path(file_path,simulated_interactions_file))[["ligand_receptor"]] 
    simulated_interactions_lst[[dataset]][[naming]] = tmp_vec 
    
    # load correct PCE file depending on the params_grid
    PCE_file = PCE_files[grepl(paste0("_" , x["FC"] , "_", ".*", x["PCE"] , ".RDS$"), PCE_files)]
    master_lst[["PCE"]][[naming]] = readRDS(file.path(file_path,PCE_file)) %>% unlist * 100 # transform to percentage
    
    
    ########################################
    ##### Generate L/R inflated per CT #####
    # Here we are iterating per dataset
    # We are also iterating across FC/PCE parameter grid
    
    # It is written in this more "messy" way for cases when we are testing for multiple sender-receiver cells
    CT_present = c("CT1","CT2")
    for(CT in CT_present)
    {
      # This creates a table with 4 columns. 
      split_table = simulated_interactions_lst[[dataset]][[naming]] %>%
        str_split_fixed("_", n = 4) %>% 
        as.data.frame()
      
      colnames(split_table) = c("ligand", "receptor1", "receptor2", "receptor3")
      
      # Find Ligand genes which were inflated in specified celltype 
      L_genes_inflated_per_CT_to_keep = split_table$ligand %>% unique %>% setdiff(.,"")
      
      # Find Receptor genes which were inflated in specified celltype inflated 
      R_genes_inflated_per_CT_to_keep = c(split_table$receptor1, split_table$receptor2, split_table$receptor3) %>% unique %>% setdiff(.,"")
      
      # Save inflated genes per CT
      if(CT == "CT1")
      {
        master_lst[[CT]][[naming]] = list(L = L_genes_inflated_per_CT_to_keep)
      } else {master_lst[[CT]][[naming]] = list(R = R_genes_inflated_per_CT_to_keep)}
      
    }
    
    ####################################################################
    ##### Generate density plots of simulated genes b/a simulation #####
    # Here we are iterating per dataset
    # We are also iterating across FC/PCE parameter grid
    
    # check if cellnames are ordered
    stopifnot(metadata$cell_ID == colnames(inflated_counts))
    
    set.seed(1)
    gene_names = master_lst[[CT]][[naming]] %>% unlist %>% sample(.,9) # sample only 9 genes as the report will contain only 9
    
    for(gene in gene_names)
    {
      # create df with counts of inflated and original count matrices
      if(gene %in% master_lst[[CT]][[naming]]$L) 
      {
        tmp_df = data.frame(inflated_counts = inflated_counts[gene,metadata$Celltype == "CT1"] %>% as.numeric,
                            original_counts = original_counts[gene,metadata$Celltype == "CT1"] %>% as.numeric) %>% melt
      } else {
        tmp_df = data.frame(inflated_counts = inflated_counts[gene,metadata$Celltype == "CT2"] %>% as.numeric,
                            original_counts = original_counts[gene,metadata$Celltype == "CT2"] %>% as.numeric) %>% melt
      }
      
      
      # calculate mean value for both counts
      mu = ddply(tmp_df, "variable", summarise, grp.mean=mean(value))
      
      p = ggplot(tmp_df, aes(x=value, color=variable)) +
        geom_density()+
        geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),
                   linetype="dashed") +
        ggtitle(paste0("density of counts | lines - mean | FC:", x["FC"], " PCE:",x["PCE"], " gene:", gene)) +
        theme_light() +
        theme(plot.title = element_text(size = 6, face = "bold"),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.title=element_blank(),
              axis.text=element_text(size=6)) 
      
      diagnostic_gene_densityPlots[[dataset]][[naming]][[CT]][[gene]] = p
    }
    
    ###########################################################
    ##### Generate Plots of real FC after semi-simulation #####
    # Here we are iterating per dataset
    # We are also iterating across FC/PCE parameter grid
    
    realFC_aftersimulation_file = file_path %>% list.files(., pattern = "realFC_aftersimulation") %>% magrittr::extract(grepl(paste0("_",x["FC"] , "_" ,".*",x["PCE"] , ".RDS$"), .))
    
    tmp_lst = readRDS(file.path(file_path,realFC_aftersimulation_file)) %>% 
      unlist %>% subset(.,!is.infinite(.)) %>% 
      plot_FCafter_semisimulation(. , theoreticalFC = x["FC"], PCE = x["PCE"])
    
    diagnostic_plots_realFC[[dataset]][[paste0("FC_",x["FC"])]][[paste0("PCE_",x["PCE"])]] = tmp_lst %>% pluck("plot")
    
    diagnostic_df_realFC[[dataset]][[paste0("index_",i)]] = tmp_lst[c(2,3,4)] %>% as.data.frame
    
    ##########################
    ##### Save AvelogCPM #####
    # filter the inflated counts based to contain cells belonging to the celltype indicated by CT_toPlot
    # This is for "compute_diagnostic_plots" function in order to save memory instead of saving entire inflated counts
    master_lst[["counts_aveLogCPM"]][[naming]] = inflated_counts[,metadata$Celltype %in% CT_present] %>% aveLogCPM
    names(master_lst[["counts_aveLogCPM"]][[naming]]) = rownames(inflated_counts)
  }
  
  ##################################
  ##### Mean var plot of genes #####
  # LR sampled are the same regarding FC and PCE parameters. Meaning that we can use any master_lst[[naming]][[CT]]
  
  genemetadata = readRDS(file.path("data/processed/", dataset,"genemetadata.RDS"))
  tmp_df_L = subset(genemetadata$mean, gene_names %in% (master_lst[[CT_present[1]]][[naming]] %>% unlist %>% as.character))
  tmp_df_R = subset(genemetadata$mean, gene_names %in% (master_lst[[CT_present[2]]][[naming]] %>% unlist %>% as.character))
  p = ggplot(genemetadata$mean, aes(x = means, y = vars)) + geom_hex() +
    geom_point(data=genemetadata$mean[genemetadata$mean$randomly_sampled,], aes(x=means, y=vars), colour="red", size=2) +
    geom_smooth(data=genemetadata$mean[genemetadata$mean$randomly_sampled,], aes(x=means, y=vars), method = "lm", color = "red") +
    xlab("variance (log10)") +
    ylab("mean (log10)") +
    geom_point(data=tmp_df_L, aes(x=means, y=vars), colour="green", size=2) +
    geom_smooth(data=tmp_df_L, aes(x=means, y=vars), method = "lm", color = "green") +
    geom_point(data=tmp_df_R, aes(x=means, y=vars), colour="yellow", size=2) +
    geom_smooth(data=tmp_df_R, aes(x=means, y=vars), method = "lm", color = "yellow") +
    ggtitle(paste0("MeanVar plot of all genes | Red - Sampled random 10 genes per bin | Green - L sampled | Yellow - R sampled | DATASET - ",dataset)) +
    theme_light()
  
  diagnostic_plots_MeanVar[[dataset]] = p
  
  # filter original avelogcpm counts to contain same genes as the inflated count matrices
  original_counts = original_counts %>% subset(rownames(original_counts) %in% rownames(inflated_counts))
  
  # generate diagnostic plots of avelogcpm and PCE 
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst, 
                                                               FC_param = c(min(FC),FC[FC <= median(FC)] %>% tail(1),max(FC)), PCE_param = c(min(PCE),PCE[PCE <= median(PCE)] %>% tail(1),max(PCE)), 
                                                               dataset = dataset, metadata = metadata, CT_toPlot = CT_present)
}

######################
##### Save plots #####
for(dataset in datasets)
{
  # save master_lst as its needed for user to see which combination of FC and PCE they should use
  saveRDS(master_lst,file.path(path_results_dir ,paste0(dataset, "_master_lst.RDS")))
  
  
  pdf(file.path(path_results_dir ,paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  
  ggarrange(plotlist = diagnostic_plots_perCT[[dataset]]$avelogcpm_fixedPCE, common.legend = TRUE) %>% print
  ggarrange(plotlist = diagnostic_plots_perCT[[dataset]]$PCE_fixedPCE) %>% print
  
  x = 1:length(diagnostic_plots_realFC[[dataset]])
  sapply(x, function(x) {ggarrange(plotlist = diagnostic_plots_realFC[[dataset]][[x]]) %>% print}) %>% print
  diagnostic_plots_MeanVar[[dataset]] %>% print
  
  
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_1_PCE_7[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_1_PCE_20[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_3_PCE_7[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_3_PCE_20[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_5_PCE_7[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_5_PCE_20[[1]], common.legend = T) %>% print
  
  dev.off()
}

##################################
##### precision/recall plots #####
##################################
statistics_results_lst_recallprecision_plot = list()
significant_interactions_lst = list()

# get all metric for all methods and datasets into a list
statistics_results_lst_recallprecision_plot = lapply(names(metric_results), function(dataset) {
  tmp_diagnostic_df_realFC = do.call(rbind, diagnostic_df_realFC[[dataset]])
  
  lapply(metric_results[[dataset]], function(method) {
    tmp_df = method %>% do.call(rbind, .) %>%
      mutate(FC_real_median = tmp_diagnostic_df_realFC$FC_real_median,
             dataset = dataset)
  }) %>% do.call(rbind.data.frame,.)
}) 

names(statistics_results_lst_recallprecision_plot) = names(metric_results)

# generate many plots
for(dataset in datasets)
{
  for(method in methods)
  {
    # Load files for upset plot
    significant_interactions_files = file.path(path_output_dir,dataset,method) %>% 
      list.files(., pattern = "significant_interactions")
    
    for(file in significant_interactions_files)
    {
      significant_interactions_lst[["all_significant"]][[method]][[file]] = read.table(file.path(path_output_dir,dataset,method,file),header = T) %>% 
        mutate(significant = as.logical(significant)) %>%
        filter(significant) %>% 
        select(ligand_receptor) %>% 
        unlist %>%
        unname
      
      significant_interactions_lst[["intersect_significant_simulated"]][[method]][[file]] = intersect(significant_interactions_lst[["all_significant"]][[method]][[file]], simulated_interactions_lst[[dataset]][[1]])
    }
  }
  
  #######################################
  ##### Generate heatmap of f1score #####
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  
  # Create a 'Combination' label for columns
  tmp_df %<>%
    mutate(Param_Comb = paste0("FC:", FC, "\nPCE:", PCE)) %>%
    select(method, Param_Comb, f1score) %>%
    pivot_wider(names_from = Param_Comb, values_from = f1score) %>%
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  # Extract metadata for annotations (the labels at the top)
  # This keeps the FC and PCE values linked to the columns
  col_meta = data.frame(colnames(tmp_df)) %>%
    separate(1, into = c("FC", "PCE"), sep = "\n") %>%
    mutate(across(everything(), ~gsub(".*:", "", .))) # Clean strings
  
  # order
  col_meta$FC = factor(col_meta$FC, levels = unique(col_meta$FC))
  col_meta$PCE = factor(col_meta$PCE, levels = unique(col_meta$PCE))
  
  # Get unique values
  unique_fcs = unique(col_meta$FC)
  unique_pces = unique(col_meta$PCE)
  
  # Using RColorBrewer for easy colors
  fc_cols = setNames(brewer.pal(length(unique_fcs), "Set1"), unique_fcs)
  pce_cols = setNames(brewer.pal(length(unique_pces), "Set2"), unique_pces)
  
  top_ann = HeatmapAnnotation(
    FC = col_meta$FC,
    PCE = col_meta$PCE,
    col = list(
      FC = fc_cols, 
      PCE = pce_cols
    )
  )
  
  # f1 score color gradient (0 to 1)
  col_f1 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  
  heatmap1 = Heatmap(tmp_df, 
                     name = "F1 Score",
                     col = col_f1,
                     top_annotation = top_ann,
                     column_split = col_meta$PCE, # Visually separate by PCE
                     cluster_rows = TRUE,     
                     cluster_columns = FALSE, 
                     width = ncol(tmp_df) * unit(6, "mm"),
                     height = nrow(tmp_df) * unit(8, "mm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 8),
                     show_column_names = FALSE,
                     rect_gp = gpar(col = "white", lwd = 1)) # Add white borders to cells
  
  #####################################################################################
  ##### Generate heatmap/upset of common significant interactions among simulated #####
  
  lst_upset_plots_FC = list()
  lst_upset_plots_PCE = list()
  lst_binarydf = list()
  for(fc in FC)
  {
    for(pce in PCE)
    {
      #tmp = lapply(significant_interactions_lst$all_significant, function(method) {method[[paste0("significant_interactions_FC_", fc, "_PercCellsExpressing_",pce,".csv")]]})
      tmp = lapply(significant_interactions_lst$intersect_significant_simulated, function(method) {method[[paste0("significant_interactions_FC_", fc, "_PercCellsExpressing_",pce,".csv")]]})
      
      # when all methods didnt retrieve a single simulated interaction
      if(length(tmp %>% unlist) == 0)
      {
        binary_df = data.frame(cellphonedb = NA,
                               singlecellsignalR = NA,
                               cellchat = NA,
                               connectome = NA,
                               natmi_specificity = NA,
                               geometric_mean = NA,
                               scseqcomm = NA,
                               log2fc = NA)
        
        lst_binarydf[[paste0("FC_" ,fc,"-PCE_",pce)]] = binary_df
      } else {
        all_strings = unique(unlist(tmp))
        
        # For each element in the list, check if the global strings exist there
        binary_df = as.data.frame(lapply(tmp, function(x) {
          as.numeric(all_strings %in% x)
        }))
        rownames(binary_df) = all_strings
        lst_binarydf[[paste0("FC_" ,fc,"-PCE_",pce)]] = binary_df
      }
      
      x = upset(
        binary_df, 
        colnames(binary_df),
        name = NULL,
        width_ratio = 0.275,
        base_annotations = list(
          'Intersection size' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = "bars") # You can style the bars here
          )
        ),
        set_sizes = (
          upset_set_size() + 
            theme(
              axis.text.x = element_text(size = 8),
              axis.title.x = element_text(size = 8)
            ) +
            # Use expand_limits to ensure the axis goes high enough for the labels
            expand_limits(y = (colSums(binary_df) %>% max) + 10) +
            geom_text(aes(label = ..count..), hjust = 1.05, stat = 'count', size = 3)
        ),
        themes = upset_default_themes(text = element_text(size = 12),
                                      legend.position = "none")
      )
      
      lst_upset_plots_FC[[paste0("FC_",fc)]][[paste0("PCE_",pce)]] = x
      lst_upset_plots_PCE[[paste0("PCE_",pce)]][[paste0("FC_",fc)]] = x
    }
  }
  
  # Function to calculate pairwise intersections for one matrix
  get_pairwise_intersects = function(df) {
    method_names = colnames(df)
    pairs = combn(method_names, 2) # Get all combinations
    
    results = apply(pairs, 2, function(p) {
      # Count rows where BOTH methods in the pair found the interaction (1)
      intersect_count = sum(df[[p[1]]] == 1 & df[[p[2]]] == 1)
      return(intersect_count)
    })
    
    names(results) = apply(pairs, 2, paste, collapse = " & ")
    return(results)
  }
  
  # Apply to all combinations of params
  summary_list = lapply(lst_binarydf, get_pairwise_intersects)
  summary_mat = do.call(cbind, summary_list)
  
  # Convert to percentage
  total_simulated = config$nLR_per_CTCTcomb # ground truth
  summary_mat_pct = (summary_mat / total_simulated) * 100
  
  # Extract metadata for annotations (the labels at the top)
  # This keeps the FC and PCE values linked to the columns
  col_meta = data.frame(colnames(summary_mat_pct)) %>%
    separate(1, into = c("FC", "PCE"), sep = "-") %>%
    mutate(FC = gsub("_",":",FC),
           PCE = gsub("_",":",PCE))
  
  # order
  col_meta$FC = factor(col_meta$FC, levels = unique(col_meta$FC))
  col_meta$PCE = factor(col_meta$PCE, levels = unique(col_meta$PCE))
  
  # Get unique values
  unique_fcs = unique(col_meta$FC)
  unique_pces = unique(col_meta$PCE)
  
  # Using RColorBrewer for easy colors
  fc_cols = setNames(brewer.pal(length(unique_fcs), "Set1"), unique_fcs)
  pce_cols = setNames(brewer.pal(length(unique_pces), "Set2"), unique_pces)
  
  top_ann = HeatmapAnnotation(
    FC = col_meta$FC,
    PCE = col_meta$PCE,
    col = list(
      FC = fc_cols, 
      PCE = pce_cols
    )
  )
  
  # f1 score color gradient (0 to 1)
  col_f1 = colorRamp2(c(0, 50, 100), c("blue", "white", "red"))
  
  heatmap2 = Heatmap(summary_mat_pct, 
                     name = "% significant int. retrieved",
                     col = col_f1,
                     cluster_rows = FALSE,
                     top_annotation = top_ann,
                     column_split = col_meta$PCE, # Visually separate by PCE
                     cluster_columns = FALSE, 
                     row_names_side = "left",
                     width = ncol(tmp_df) * unit(6, "mm"),
                     height = nrow(tmp_df) * unit(20, "mm"),
                     column_names_gp = gpar(fontsize = 8),, 
                     row_names_gp = gpar(fontsize = 10), 
                     show_column_names = FALSE,
                     rect_gp = gpar(col = "white", lwd = 1)) # Add white borders to cells
  
  ###########################################
  ##### Generate alternative to heatmap #####
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  
  # Calculate the mean F1 score per method to order them logically
  method_order = tmp_df %>% 
    group_by(method) %>% 
    dplyr::summarise(m = mean(f1score)) %>% 
    arrange(m) %>% 
    pull(method)
  
  tmp_df$method = factor(tmp_df$method, levels = method_order)
  
  tmp_df$PCE_label = factor(tmp_df$PCE, levels = c(4, 7, 10, 20, 30, 40))
  
  heatmap1_alternative = ggplot(tmp_df, aes(x = f1score, y = method)) +
    # Add a light background line for each method
    geom_segment(aes(x = 0, xend = 1, y = method, yend = method), 
                 color = "gray90", size = 0.5) +
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
    labs(
      title = "Comparing F1 scores across sparsity (PCE) and signal (FC)",
      x = "F1 Score",
      y = "",
      color = "Fold Change"
    )
  
  ######################
  ##### UpSet plot #####
  
  ##### all_significant_interactions
  all_significant_interactions = lapply(significant_interactions_lst[["all_significant"]] , function(method) {
    # find all unique elements
    all_strings = unique(unlist(significant_interactions_lst[["all_significant"]]))
    # For each element in the list, check if the global strings exist there
    binary_df = as.data.frame(lapply(method, function(x) {
      as.numeric(all_strings %in% x)
    }))
    
    data.frame(method = ifelse(rowSums(binary_df) > 0, 1,0), row.names = all_strings)
  }) %>% do.call(cbind,.)
  
  colnames(all_significant_interactions) = names(significant_interactions_lst[["all_significant"]])
  
  upset1 = upset(
    all_significant_interactions, 
    colnames(all_significant_interactions),
    sort_intersections_by = 'degree',
    sort_intersections = 'descending', # 'descending' puts single-method hits on the right
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        mapping = aes(fill = "bars") # You can style the bars here
      )
    ),
    set_sizes = (
      upset_set_size() + 
        # Use expand_limits to ensure the axis goes high enough for the labels
        expand_limits(y = (colSums(all_significant_interactions) %>% max) + 10) +
        geom_text(aes(label = ..count..), hjust = 1.2, stat = 'count', size = 3)
    ),
    themes = upset_default_themes(text = element_text(size = 12),
                                  legend.position = "none")
  )
  
  ##### intersect_significant_simulated_interactions
  intersect_significant_simulated_interactions = lapply(significant_interactions_lst[["intersect_significant_simulated"]] , function(method) {
    # find all unique elements
    all_strings = simulated_interactions_lst[[dataset]][[1]]
    # For each element in the list, check if the global strings exist there
    binary_df = as.data.frame(lapply(method, function(x) {
      as.numeric(all_strings %in% x)
    }))
    
    data.frame(method = ifelse(rowSums(binary_df) > 0, 1,0), row.names = all_strings)
  }) %>% do.call(cbind,.)
  
  colnames(intersect_significant_simulated_interactions) = names(significant_interactions_lst[["intersect_significant_simulated"]])
  
  upset2 = upset(
    intersect_significant_simulated_interactions, 
    colnames(intersect_significant_simulated_interactions),
    sort_intersections_by = 'degree',
    sort_intersections = 'descending', # 'descending' puts single-method hits on the right
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        mapping = aes(fill = "bars") # You can style the bars here
      )
    ),
    set_sizes = (
      upset_set_size() + 
        # Use expand_limits to ensure the axis goes high enough for the labels
        expand_limits(y = (colSums(intersect_significant_simulated_interactions) %>% max) + 10) +
        geom_text(aes(label = ..count..), hjust = 1.1, stat = 'count')
    ),
    themes = upset_default_themes(text = element_text(size = 12),
                                  legend.position = "none")
  )
  
  ##########################################
  ##### Generate Precision recall plot #####
  
  tmp_PCE = statistics_results_lst_recallprecision_plot[[dataset]]$PCE
  tmp_FC = statistics_results_lst_recallprecision_plot[[dataset]]$FC
  
  # Plot precision recall curves when generating 1 plot per PCE
  lst_precision_recall_byPCE = list()
  for(pce in tmp_PCE)
  {
    data = filter(statistics_results_lst_recallprecision_plot[[dataset]], PCE == pce)
    p = ggplot(data, aes(y = precision, x = recall , color = method)) + 
      geom_point(size = 1.5) + 
      geom_line() +
      theme_light() +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))  + 
      theme(legend.title = element_text(size = 7), 
            legend.text = element_text(size = 7),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),  
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            plot.title = element_text(size=10)) +
      ggtitle(paste0(paste0("PCE_",pce)))
    lst_precision_recall_byPCE[[paste0("PCE_",pce)]] = p
  }
  
  # Plot precision recall curves when generating 1 plot per PCE
  lst_precision_recall_byFC = list()
  for(fc in tmp_FC)
  {
    data = filter(statistics_results_lst_recallprecision_plot[[dataset]], FC == fc)
    p = ggplot(data, aes(y = precision, x = recall , color = method)) + 
      geom_point(size = 1.5) + 
      geom_line() + 
      theme_light() +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))  + 
      theme(legend.title = element_text(size = 7), 
            legend.text = element_text(size = 7),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),  
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            plot.title = element_text(size=10)) +
      ggtitle(paste0(paste0("FC_",fc)))
    lst_precision_recall_byFC[[paste0("FC_",fc)]] = p
  }
  
  ##################################
  ##### Generate f1 score plot #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  p1 = ggplot(tmp_df, aes(x = PCE, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~FC ) +
    ggtitle("Faceted by FC") +
    theme_light()
  
  p2 = ggplot(tmp_df, aes(x = FC, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~PCE) +
    ggtitle("Faceted by PCE") +
    theme_light()
  
  ##########################################
  ##### Generate ranking LR genes plot #####
  p3 = ggplot(tmp_df, aes(x = FC, y = normalizedRank_simulated_among_significant, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~PCE) +
    ggtitle("Faceted by PCE") +
    theme_light()
  
  ##### Save plots
  pdf(file.path(path_results_dir ,paste0(dataset, "_results_plots.pdf")), width = 15, height = 10)
  wrap_plots(lst_precision_recall_byPCE, guides = "collect")  %>% print
  wrap_plots(lst_precision_recall_byFC, guides = "collect") %>% print
  p1 %>% print
  p2 %>% print
  p3 %>% print
  heatmap1 %>% print
  heatmap1_alternative %>% print
  heatmap2 %>% print
  (wrap_plots(lst_upset_plots_FC$FC_0.3) + 
    plot_annotation( title = paste0("FC = 1 and PCE = ", paste(PCE, collapse = ", ")))) %>% print
  (wrap_plots(lst_upset_plots_FC$FC_3) + 
    plot_annotation(title = paste0("FC = 5 and PCE = ", paste(PCE, collapse = ", ")))) %>% print
  (wrap_plots(lst_upset_plots_PCE$PCE_7) + 
    plot_annotation(title = paste0("PCE = 10 and FC = ", paste(FC, collapse = ", ")))) %>% print
  (wrap_plots(lst_upset_plots_PCE$PCE_20) + 
    plot_annotation(title = paste0("PCE = 20 and FC = ", paste(FC, collapse = ", ")))) %>% print
  
  upset1 %>% print
  grid.text("All significant interactions retrieved across any combination of params", 
            x = 0.2, y = 0.98, 
            gp = gpar(fontsize = 8, fontface = "bold"))
  
  upset2 %>% print
  grid.text("Common significant/simulated interactions retrieved across any combination of params", 
            x = 0.2, y = 0.98, 
            gp = gpar(fontsize = 8, fontface = "bold"))
  
  dev.off()
}

###################################
##### f1score across datasets #####

x = do.call(rbind.data.frame,statistics_results_lst_recallprecision_plot) %>%
  mutate(
    PCE = as.numeric(as.character(PCE)),
    FC = as.numeric(as.character(FC))
  ) %>%
  # 2. Sort the actual data rows by these numbers
  arrange(PCE, FC) %>%
  # 3. Create a unique parameter string that preserves this order
  mutate(params = paste0("PCE_", PCE, "_FC_", FC)) %>%
  # 4. CRITICAL: Convert to factor using the sorted unique values
  mutate(params = factor(params, levels = unique(params)))

# 1. Calculate relative scores within each group
df_relative <- x %>%
  group_by(dataset, params) %>%
  dplyr::mutate(relative_score = f1score / sum(f1score, na.rm = TRUE)) %>%
  dplyr::summarise(
    # Store the relative scores as a named vector
    rel_vec = list(setNames(relative_score, method)),
    .groups = "drop"
  )

# 2. Pivot to the Matrix of Lists
mat_rel <- df_relative %>%
  pivot_wider(names_from = params, values_from = rel_vec) %>%
  tibble::column_to_rownames("dataset") %>%
  as.matrix()

# 3. Method Colors (Ensure you have 8 distinct colors)
method_cols <- setNames(brewer.pal(8, "Set2"), unique(x$method))

# 1. Extract the PCE values from the column names of your matrix
# If your column names are "PCE_4_FC_5", this grabs the "4"
column_pce_values <- str_split_i(colnames(mat_rel), "_", 2)

# 2. (Optional) Convert to factor to ensure they stay in numerical order (4, 7, 10...)
column_pce_values <- factor(column_pce_values, levels = unique(column_pce_values))


ht <- Heatmap(mat_rel, 
              show_heatmap_legend = FALSE,
              col = method_cols,
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              column_title = "Relative f1score",
              # --- THE SPLIT LOGIC ---
              column_split = column_pce_values, 
              column_gap = unit(5, "mm"),        # Adjust the width of the gap here
              border = TRUE,                     # Adds a border around each split block
              rect_gp = gpar(type = "none",         # Keeps your custom bars
                             col = "black",         # Border color
                             lwd = 2),              # INCREASE THIS (e.g., 2 or 3)
              
              cell_fun = function(j, i, x, y, width, height, fill) {
                scores <- mat_rel[[i, j]]
                if(is.null(scores) || any(is.na(scores))) return()
                
                # 'current_y' tracks the bottom of the next segment in the stack
                current_y = y - height/2
                
                for(k in seq_along(scores)) {
                  seg_h = height * scores[k] # height is proportional to relative_score
                  
                  grid.rect(
                    x = x, 
                    y = current_y + seg_h/2, 
                    width = width * 0.9, # 0.9 adds a small gap between columns
                    height = seg_h,
                    gp = gpar(fill = method_cols[names(scores)[k]], col = NA)
                  )
                  # Move the starting point up for the next method
                  current_y = current_y + seg_h
                }
              })

pdf(file.path(path_results_dir ,"results_acrossdatasets_plots.pdf"), width = 12, height = 10)
# Draw with Legend
draw(ht, annotation_legend_list = list(Legend(labels = names(method_cols), 
                                              title = "Methods", 
                                              legend_gp = gpar(fill = method_cols))))

dev.off()

####################################################
##### Create recommendation hierarchical tree #####
# it was done manually.

"""
# select only 3 FC and 3 PCE for the hierarchical plot
FC = c(0.3,1,10)
PCE = c(7,20,40)

library(scater) 
out = list()
for(dataset in datasets)
{
  # Change the file according to your dataset
  master_lst = readRDS(file.path(path_results_dir ,paste0(dataset, '_master_lst.RDS')))
  benchmark_simulated_genes = c(master_lst$CT1[[1]]$L , master_lst$CT2[[1]]$R) %>% unique # simulated genes are the same for entire dataset
  
  out[[dataset]] = compare_expressionProfile(dataset = dataset,
                                  FC_values = FC,
                                  PCE_values = PCE,
                                  benchmark_simulated_genes = benchmark_simulated_genes)
  
  #pdf(file.path(path_results_dir ,paste0(dataset,'_method_recommendation_plots.pdf')), width = 12, height = 10)
  #ggarrange(plotlist = out$plots$PCE7, common.legend = T) %>% print
  #ggarrange(plotlist = out$plots$PCE20, common.legend = T) %>% print
  #ggarrange(plotlist = out$plots$PCE40, common.legend = T) %>% print
  #dev.off()
}

df_10x = lapply(out$`10x`$quantiles, function(x) {
  as.data.frame(x)['50%',]
}) %>%
  bind_rows(.id = 'PCE') %>%
  pivot_longer(
    cols = starts_with('FC'), 
    names_to = 'FC', 
    values_to = 'Expression'
  ) %>%
  mutate(
    technology = dataset,
    PCE = str_remove(PCE, 'PCE'),
    FC = str_remove(FC, 'FC')
  )

df_smartseq2 = lapply(out$SMARTseq2$quantiles, function(x) {
  as.data.frame(x)['50%',]
}) %>%
  bind_rows(.id = 'PCE') %>%
  pivot_longer(
    cols = starts_with('FC'), 
    names_to = 'FC', 
    values_to = 'Expression'
  ) %>%
  mutate(
    technology = dataset,
    PCE = str_remove(PCE, 'PCE'),
    FC = str_remove(FC, 'FC')
  )

df_vasaseq = lapply(out$VASAseq$quantiles, function(x) {
  as.data.frame(x)['50%',]
}) %>%
  bind_rows(.id = 'PCE') %>%
  pivot_longer(
    cols = starts_with('FC'), 
    names_to = 'FC', 
    values_to = 'Expression'
  ) %>%
  mutate(
    technology = dataset,
    PCE = str_remove(PCE, 'PCE'),
    FC = str_remove(FC, 'FC')
  )

df_10x_immune = lapply(out$`10x_immune_R1`$quantiles, function(x) {
  as.data.frame(x)['50%',]
}) %>%
  bind_rows(.id = 'PCE') %>%
  pivot_longer(
    cols = starts_with('FC'), 
    names_to = 'FC', 
    values_to = 'Expression'
  ) %>%
  mutate(
    technology = dataset,
    PCE = str_remove(PCE, 'PCE'), 
    FC = str_remove(FC, 'FC')
  )

df_BD = lapply(out$BD_wholeTranscriptome_P1$quantiles, function(x) {
  as.data.frame(x)['50%',]
}) %>%
  bind_rows(.id = 'PCE') %>%
  pivot_longer(
    cols = starts_with('FC'), 
    names_to = 'FC', 
    values_to = 'Expression'
  ) %>%
  mutate(
    technology = dataset,
    PCE = str_remove(PCE, 'PCE'),
    FC = str_remove(FC, 'FC')
  )


# check the f1scores
tmp_df = statistics_results_lst_recallprecision_plot$`10x_immune_R1`

# Create a 'Combination' label for columns
tmp_df %<>%
  mutate(Param_Comb = paste0('FC:', FC, '\nPCE:', PCE)) %>%
  select(method, Param_Comb, f1score) %>%
  pivot_wider(names_from = Param_Comb, values_from = f1score) %>%
  tibble::column_to_rownames('method') %>%
  as.matrix()

x = tmp_df[,grepl('PCE:40', colnames(tmp_df))]
x[,1:3] %>% rowMeans() %>% sort(decreasing = T)
x[,4:6] %>% rowMeans() %>% sort(decreasing = T)

############## 
df_smartseq2 = data.frame(technology = 'SMARTseq2',
                          PCE = c(7,7,20,20,40,40),
                          Expression = c('<4.45','>4.45','<3.87','>3.87','<3.62','>3.62'),
                          Method = c('log2fc','log2fc',
                                     'log2fc|geometric_mean|natmi_specificity|cellphonedb','log2fc|geometric_mean|cellphonedb|natmi_specificity',
                                     'cellchat','singlecellsignalR|cellchat'))

df_10x = data.frame(technology = '10x_1',
                    PCE = c(7,7,20,20,40,40),
                    Expression = c('<1.6','>1.6','<1.42','>1.42','<1.34','>1.34'),
                    Method = c('geometric_mean|cellphonedb','geometric_mean|cellphonedb',
                               'scseqcomm|singlecellsignalR','singlecellsignalR',
                               'scseqcomm|singlecellsignalR|cellchat','scseqcomm|singlecellsignalR|cellchat'))

df_vasaseq = data.frame(technology = 'VASAseq',
                        PCE = c(7,7,20,20,40,40),
                        Expression = c('<2.42','>2.42','<2.38','>2.38','<2.45','>2.45'),
                        Method = c('-','log2fc',
                                   '-','log2fc',
                                   'singlecellsignalR','log2fc|singlecellsignalR'))

df_10x_immune = data.frame(technology = '10x_2',
                        PCE = c(7,7,20,20,40,40),
                        Expression = c('<1.22','>1.22','<2.38','>2.38','<2.45','>2.45'),
                        Method = c('geometric_mean|cellphonedb','geometric_mean|cellphonedb|singlecellsignalR',
                                   'scseqcomm|singlecellsignalR','scseqcomm|singlecellsignalR',
                                   'scseqcomm|singlecellsignalR|cellchat','scseqcomm|singlecellsignalR|cellchat'))

df_BD = data.frame(technology = 'BDRhapsody',
                        PCE = c(7,7,20,20,40,40),
                        Expression = c('<2.42','>2.42','<2.38','>2.38','<2.45','>2.45'),
                        Method = c('geometric_mean|cellphonedb','geometric_mean|cellphonedb',
                                   'singlecellsignalR|geometric_mean|cellphonedb','singlecellsignalR',
                                   'singlecellsignalR|cellchat','singlecellsignalR|cellchat'))

df = bind_rows(list(df_smartseq2,df_10x,df_vasaseq,df_10x_immune,df_BD))

# create by gemini
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

library(DiagrammeR)

library(DiagrammeR)

graph <- grViz('
digraph selection_tree_updated {
  # Global Settings
  graph [layout = dot, ranksep = 0.5, nodesep = 0.1]
  node [shape = box, fontname = Helvetica, style = filled, fillcolor = '#FDFEFE']
  edge [color = '#566573', arrowhead = vee, fontsize = 10]
  
  # ROOT
  Start [label = 'Sequencing Technology', fillcolor = '#D6EAF8', style = 'filled,bold']
  
  # TECHNOLOGY NODES
  node [fillcolor = '#AED6F1']
  T_SMART [label = 'SMARTseq2']
  T_10x1  [label = '10x_1']
  T_10x2  [label = '10x_2']
  T_VASA  [label = 'VASAseq']
  T_BD    [label = 'BDRhapsody']
  
  # PCE NODES
  node [fillcolor = '#FCF3CF', height = 0.3, width = 0.6]
  S7 [label='PCE 7']; S20 [label='PCE 20']; S40 [label='PCE 40']
  X1_7 [label='PCE 7']; X1_20 [label='PCE 20']; X1_40 [label='PCE 40']
  X2_7 [label='PCE 7']; X2_20 [label='PCE 20']; X2_40 [label='PCE 40']
  V7 [label='PCE 7']; V20 [label='PCE 20']; V40 [label='PCE 40']
  B7 [label='PCE 7']; B20 [label='PCE 20']; B40 [label='PCE 40']
  
  # METHOD NODES (Color coded by Technology logic)
  node [fillcolor = '#D5F5E3', fontsize = 9]
  M_log2 [label = 'log2fc']
  M_S20  [label = 'log2fc\ngeometric_mean\nnatmi_specificity\ncellphonedb']
  M_chat [label = 'cellchat']
  M_S40  [label = 'singlecellsignalR\ncellchat']
  
  node [fillcolor = '#E8DAEF']
  M_G_C  [label = 'geometric_mean\ncellphonedb']
  M_X20  [label = 'scseqcomm\nsinglecellsignalR']
  M_sigR [label = 'singlecellsignalR']
  M_X40  [label = 'scseqcomm\nsinglecellsignalR\ncellchat']
  
  node [fillcolor = '#FADBD8']
  M_None [label = 'No recommendation (-)']
  
  node [fillcolor = '#FEF9E7']
  M_I7   [label = 'geometric_mean\ncellphonedb\nsinglecellsignalR']
  M_B20  [label = 'singlecellsignalR\ngeometric_mean\ncellphonedb']
  M_B40  [label = 'singlecellsignalR\ncellchat']
  M_V40  [label = 'log2fc\nsinglecellsignalR']
  
  # CONNECTIONS
  Start -> {T_SMART T_10x1 T_10x2 T_VASA T_BD}
  
  # SMARTseq2 Logic
  T_SMART -> {S7 S20 S40}
  S7  -> M_log2 [label = 'Any']
  S20 -> M_S20  [label = 'Any']
  S40 -> M_chat [label = '< 3.62']
  S40 -> M_S40  [label = '> 3.62']
  
  # 10x_1 Logic (Formerly 10x Genomics)
  T_10x1 -> {X1_7 X1_20 X1_40}
  X1_7  -> M_G_C  [label = 'Any']
  X1_20 -> M_X20  [label = '< 1.42']
  X1_20 -> M_sigR [label = '> 1.42']
  X1_40 -> M_X40  [label = 'Any']
  
  # 10x_2 Logic (Formerly 10x Immune)
  T_10x2 -> {X2_7 X2_20 X2_40}
  X2_7  -> M_G_C  [label = '< 1.22']
  X2_7  -> M_I7   [label = '> 1.22']
  X2_20 -> M_X20  [label = 'Any']
  X2_40 -> M_X40  [label = 'Any']
  
  # VASAseq Logic
  T_VASA -> {V7 V20 V40}
  V7  -> M_None [label = '< 2.42']
  V7  -> M_log2 [label = '> 2.42']
  V20 -> M_None [label = '< 2.38']
  V20 -> M_log2 [label = '> 2.38']
  V40 -> M_sigR [label = '< 2.45']
  V40 -> M_V40  [label = '> 2.45']
  
  # BDRhapsody Logic
  T_BD -> {B7 B20 B40}
  B7  -> M_G_C  [label = 'Any']
  B20 -> M_B20  [label = '< 2.38']
  B20 -> M_sigR [label = '> 2.38']
  B40 -> M_B40  [label = 'Any']
}
')

graph

svg_text = export_svg(graph)
rsvg_pdf(charToRaw(svg_text), file = file.path(path_results_dir ,'reccomendation_hierarchincalTree.pdf'))

"""