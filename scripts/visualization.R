# Load package
library(dplyr)
library(optparse)
library(ggplot2)
library(jsonlite)
library(ggpubr)
library(stringr)
library(edgeR)
library(ggrepel)
library(reshape2)
library(purrr) 
library(magrittr)
library(plyr)
library(tidyr)
library(ComplexUpset)
library(ComplexHeatmap)
source("scripts/helper_functions.R")

# This function was entirely taken from package UpSetR
fromList = function (input) 
{
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  return(data)
}

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

config = yaml::read_yaml(path_config.yaml)

# Setting parameters
methods = config$methods %>% unlist
datasets = config$datasets %>% unlist
FC = config$semiSimulation$FC %>% unlist %>% as.double
PCE = config$semiSimulation$perc_cells_expressing %>% unlist %>% as.integer

############################
##### Diagnostic plots #####
############################

diagnostic_plots_lst = list()
diagnostic_plots_perCT= list()
diagnostic_plots_MeanVar = list()
diagnostic_plots_realFC = list()
diagnostic_df_realFC = list()
diagnostic_gene_densityPlots = list()
UpSet_plot_lst = list()
for(dataset in datasets)
{
  # load files
  file_path = file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB"))
  inflated_counts_files = file_path %>% list.files(., pattern = "sc_inflated_counts")
  original_counts = read.table(file.path("data/",dataset, "raw_counts.tsv"))  # load original counts
  rownames(original_counts) = rownames(original_counts) %>% toupper()
  colnames(original_counts) = gsub("[.-]","_" , colnames(original_counts))
  
  PCE_files = file_path %>% list.files(., pattern = "perc_cells_expressing_perGene")
  simulated_interactions_files = file_path %>% list.files(., pattern = "simulated_interactions")
  
  # list to save all data
  master_lst_diagnosticPlots = list()
  
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
    tmp_vec = readRDS(file.path(file_path,simulated_interactions_file))[[1]] %>% magrittr::extract(!str_detect(. , "subunit")) # remove the subunit genes that we inflated -> if we dont remove them, we will have much more FN
    master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] = tmp_vec 
    rm(tmp_vec)
    
    # load correct PCE file depending on the params_grid
    PCE_file = PCE_files[grepl(paste0("_" , x["FC"] , "_", ".*", x["PCE"] , ".RDS$"), PCE_files)]
    master_lst_diagnosticPlots[[naming]][["PCE"]] = readRDS(file.path(file_path,PCE_file)) %>% unlist * 100 # transform to percentage
    
    #################################################################################
    ##### Generate UpSet plot comparing significant interactions across methods #####
    # Here we are iterating per dataset
    # We are also iterating across FC/PCE parameter grid
    
    # Iterate over methods to get the all the True Positives that methods picked up
    for(method in methods)
    {
      significant_interactions_files = file.path(path_output_dir,dataset,method) %>%
        list.files(., pattern = "significant_interactions")
      significant_interactions_file = significant_interactions_files %>% magrittr::extract(grepl(paste0("_" , x["FC"] , "_",".*",x["PCE"] , ".RDS$"), significant_interactions_files))
      
      tmp_vec = readRDS(file.path(path_output_dir,dataset,method,significant_interactions_file))$CT1_CT2 %>% 
        mutate(ligand_receptor = str_c(.$ligand ,"_", .$receptor)) %>%
        select(ligand_receptor) %>%
        unlist %>%
        unname 
      
      tmp_vec2 = intersect(tmp_vec, master_lst_diagnosticPlots[[naming]][["simulated_interactions"]])
      tmp_vec3 = setdiff(tmp_vec, master_lst_diagnosticPlots[[naming]][["simulated_interactions"]])
      
      master_lst_diagnosticPlots[[naming]][["TP_byMethod"]][[method]] = tmp_vec2
      master_lst_diagnosticPlots[[naming]][["FP_byMethod"]][[method]] = tmp_vec3
    }
    
    # generate upset plot 
    # check if there are not TP
    if(!identical(tmp_vec2, character(0)))
    {
      UpSet_plot_lst[[dataset]][[naming]][["UpSet_plot_TP_byMethod"]] = ComplexUpset::upset(fromList(master_lst_diagnosticPlots[[naming]][["TP_byMethod"]]), intersect = methods) + ggtitle(paste("Intersection of TP -",naming))
    } else {UpSet_plot_lst[[dataset]][[naming]][["UpSet_plot_TP_byMethod"]] = ggplot(data.frame()) + ggtitle("No TP for this method") + theme_bw()}
    UpSet_plot_lst[[dataset]][[naming]][["UpSet_plot_FP_byMethod"]] = ComplexUpset::upset(fromList(master_lst_diagnosticPlots[[naming]][["FP_byMethod"]]), intersect = methods) + ggtitle(paste("Intersection of FP -",naming))
    rm(tmp_vec, tmp_vec2, tmp_vec3)
    ########################################
    ##### Generate L/R inflated per CT #####
    # Here we are iterating per dataset
    # We are also iterating across FC/PCE parameter grid
    
    # It is written in this more "messy" way for cases when we are testing for multiple sender-receiver cells
    CT_present = c("CT1","CT2")
    for(CT in CT_present)
    {
      # Find Ligand genes which were inflated in specified celltype 
      L_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] %>% str_split(.,"_") %>% 
        lapply("[[", 1) %>% unlist
      
      # Find Receptor genes which were inflated in specified celltype inflated 
      R_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] %>% str_split(.,"_") %>% 
        lapply("[[", 2) %>% unlist
      
      # Save inflated genes per CT
      master_lst_diagnosticPlots[[naming]][[CT]] = list(L = L_genes_inflated_per_CT_to_keep %>% unique , R = R_genes_inflated_per_CT_to_keep %>% unique)
    }
    
    ####################################################################
    ##### Generate density plots of simulated genes b/a simulation #####
    # Here we are iterating per dataset
    # We are also iterating across FC/PCE parameter grid
    
    # check if cellnames are ordered
    stopifnot(metadata$cell_ID == colnames(inflated_counts))
    
    set.seed(1)
    gene_names = master_lst_diagnosticPlots[[naming]][[CT]] %>% unlist %>% sample(.,12) # sample only 12 genes as the report will contain only 12
    
    for(gene in gene_names)
    {
      # create df with counts of inflated and original count matrices
      if(gene %in% master_lst_diagnosticPlots[[naming]][[CT]]$L) 
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
    master_lst_diagnosticPlots[[naming]][["counts_aveLogCPM"]] = inflated_counts[,metadata$Celltype %in% CT_present[1]] %>% aveLogCPM
    names(master_lst_diagnosticPlots[[naming]][["counts_aveLogCPM"]]) = rownames(inflated_counts)
  }
    
  ##################################
  ##### Mean var plot of genes #####
  # LR sampled are the same regarding FC and PCE parameters. Meaning that we can use any master_lst_diagnosticPlots[[naming]][[CT]]
  
  genemetadata = readRDS(file.path("data/processed/", dataset,"genemetadata.RDS"))
  tmp_df_L = subset(genemetadata$mean, gene_names %in% (master_lst_diagnosticPlots[[naming]][[CT_present[1]]] %>% unlist %>% as.character))
  tmp_df_R = subset(genemetadata$mean, gene_names %in% (master_lst_diagnosticPlots[[naming]][[CT_present[2]]] %>% unlist %>% as.character))
  p = ggplot(genemetadata$mean, aes(x = means, y = vars)) + geom_hex() +
    geom_point(data=genemetadata$mean[genemetadata$mean$gene_to_use,], aes(x=means, y=vars), colour="red", size=2) +
    geom_smooth(data=genemetadata$mean[genemetadata$mean$gene_to_use,], aes(x=means, y=vars), method = "lm", color = "red") +
    xlab("variance (log10)") +
    ylab("mean (log10)") +
    geom_point(data=tmp_df_L, aes(x=means, y=vars), colour="green", size=2) +
    geom_smooth(data=tmp_df_L, aes(x=means, y=vars), method = "lm", color = "green") +
    geom_point(data=tmp_df_R, aes(x=means, y=vars), colour="yellow", size=2) +
    geom_smooth(data=tmp_df_R, aes(x=means, y=vars), method = "lm", color = "yellow") +
    ggtitle(paste0("MeanVar plot of all genes | Red - Sampled random 10 genes per bin | Green - L sampled | Yellow - R sampled | DATASET - ",dataset))
  
  diagnostic_plots_MeanVar[[dataset]] = p
  
  # filter original avelogcpm counts to contain same genes as the inflated count matrices
  original_counts = original_counts %>% subset(rownames(original_counts) %in% rownames(inflated_counts))
  
  # generate diagnostic plots of avelogcpm and PCE 
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
                                                               FC_param = c(min(FC),median(FC),max(FC)), PCE_param = c(min(PCE),median(PCE),max(PCE)), dataset = dataset, metadata = metadata, CT_toPlot = CT_present[1])
}

######################
##### Save plots #####
for(dataset in datasets)
{
  pdf(file.path(path_results_dir ,paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  #do.call(ggarrange,c(diagnostic_plots_lst[[dataset]]$avelogcpm_fixedPCE, common.legend = TRUE)) %>% print
  #do.call(ggarrange,diagnostic_plots_lst[[dataset]]$PCE_fixedPCE) %>% print
  ggarrange(plotlist = diagnostic_plots_perCT[[dataset]]$avelogcpm_fixedPCE, common.legend = TRUE) %>% print
  ggarrange(plotlist = diagnostic_plots_perCT[[dataset]]$PCE_fixedPCE) %>% print
  
  x = 1:length(diagnostic_plots_realFC[[dataset]])
  sapply(x, function(x) {ggarrange(plotlist = diagnostic_plots_realFC[[dataset]][[x]]) %>% print}) %>% print
  diagnostic_plots_MeanVar[[dataset]] %>% print
  
  
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_1_PCE_10[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_1_PCE_50[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_3_PCE_10[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_3_PCE_50[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_5_PCE_10[[1]], common.legend = T) %>% print
  ggarrange(plotlist = diagnostic_gene_densityPlots[[dataset]]$FC_5_PCE_50[[1]], common.legend = T) %>% print
  
  dev.off()
}

##################################
##### precision/recall plots #####
##################################

# plot TPR/sensitivity/recall
statistics_results_lst_recallprecision_plot = list()
master_lst_precision_recall = list()
ranking_LRgenes_lst_plot = list()
for(dataset in datasets)
{
  for(method in methods)
  {
    CT_statistics_files = file.path(path_output_dir,dataset,"metrics/",method) %>%
      list.files(., pattern = "CT_statistics")
    master_lst_precision_recall = list()
    for(file in CT_statistics_files)
    {
      master_lst_precision_recall[[file]] = read.csv((file.path(path_output_dir,dataset,"metrics/",method,file))) %>% unlist
    }
    
    #############################################################################
    ##### Ranking of LR genes across top n of significant hits from methods #####
    
    tmp_ranking_LRgenes_lst = list()
    ranking_LRgenes_files = file.path(path_output_dir,dataset,"metrics/",method) %>% list.files(., pattern = "ranking_LRgenes")
    
    for(file in ranking_LRgenes_files)
    {
      tmp_ranking_LRgenes_lst[[file]] = read.csv((file.path(path_output_dir,dataset,"metrics/",method,file))) 
    }
    ranking_LRgenes_lst_plot[[dataset]][[method]] = do.call(rbind,tmp_ranking_LRgenes_lst) %>% mutate(. , method = method)
    
    rm(tmp_ranking_LRgenes_lst)
    
    ##### replace the theoretical FC for the median of all effective FC for all genes
    tmp_diagnostic_df_realFC = do.call(rbind, diagnostic_df_realFC[[dataset]])
    
    # Add FC and PCE columns
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] = do.call(rbind,master_lst_precision_recall) %>% as.data.frame
    
    x = rownames(statistics_results_lst_recallprecision_plot[[dataset]][[method]]) %>% stringr::str_split("_")
    tmp_FC = lapply(x, "[[" , 4) %>% unlist %>% as.double
    PCE = lapply(x, "[[" , 6) %>% unlist %>% as.double
    
    # replace the theoretical FC for the median of effective FC across genes
    FC_real = vector()
    for(i in 1:length(PCE))
    {
      FC_real = append(FC_real,filter(tmp_diagnostic_df_realFC, theoreticalFC == tmp_FC[i] & PCE_real == PCE[i]) %>% 
                         select("FC_real_median") %>% unlist %>% round(1)) 
    }
    rm(tmp_diagnostic_df_realFC)
    
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] %<>% mutate(.,FC = FC_real, PCE = PCE, method = method)
  }

  # for precision recall plots
  statistics_results_lst_recallprecision_plot[[dataset]] = do.call(rbind, statistics_results_lst_recallprecision_plot[[dataset]]) %>% mutate(., dataset = dataset)
  
  # Add FC and PCE info for ranking LR genes plot 
  ranking_LRgenes_lst_plot[[dataset]] = ranking_LRgenes_lst_plot[[dataset]] %>% do.call(rbind,.) %>%
    mutate(PCE = statistics_results_lst_recallprecision_plot[[dataset]]$PCE ,
           FC = statistics_results_lst_recallprecision_plot[[dataset]]$FC)
  #######################################
  ##### Generate heatmap of f1score #####
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  tmp_df$Parameters = paste0("FC_",tmp_df$FC,"_PCE_",tmp_df$PCE) # create additional column with parameter combination
  heatmap1 = ggplot(tmp_df, aes(method, Parameters, fill= f1score)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("heatmap of f1score sorted by increasing values of FC")
  
  tmp_df$Parameters = paste0("PCE_",tmp_df$PCE,"_FC_",tmp_df$FC) # create additional column with parameter combination
  heatmap2 = ggplot(tmp_df, aes(method, Parameters, fill= f1score)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("heatmap of f1score sorted by increasing values of PCE")
  
  # create xtabs table
  heatmap_mtx = xtabs(f1score ~ ., tmp_df %>% select(c(Parameters, f1score, method)))
  x = rownames(heatmap_mtx) %>% str_split("_")
  FC = x %>% lapply(.,"[[",2) %>% unlist %>% as.numeric
  PCE = x %>% lapply(.,"[[",4) %>% unlist %>% as.numeric
  
  # generate annotations and heatmap itself
  ha = rowAnnotation(PCE = PCE)
  ha2 = rowAnnotation(FC = FC,
                      col = setNames(RColorBrewer::brewer.pal(name = "Spectral", n = 9), FC) %>% as.list)
  library(circlize)
  col_fun = colorRamp2(c(min(heatmap_mtx), max(heatmap_mtx)), c("white", "coral"))
  
  Heatmap(heatmap_mtx, 
          left_annotation = rowAnnotation(PCE = anno_block(gp = gpar(fill = 2:4),
                                                           labels = paste0("PCE_",unique(PCE)) ,
                                                           labels_gp = gpar(col = "white", fontsize = 10))), 
          right_annotation = ha2,
          split = PCE,
          col = col_fun,
          show_row_dend = F,
          show_column_dend = F, 
          show_row_names = F,
          name = "f1score",
          column_names_rot = 60,
          column_title = "f1score of CCC methods extracting added signal")
  
  rm(tmp_df)
  ##########################################
  ##### Generate Precision recall plot #####
  
  PCE = as.factor(statistics_results_lst_recallprecision_plot[[dataset]]$PCE)
  FC = as.factor(statistics_results_lst_recallprecision_plot[[dataset]]$FC)

  # PLot precision recall curves when generating 1 plot per PCE
  lst_precision_recall_byPCE = list()
  for(pce in levels(PCE))
  {
    data = filter(statistics_results_lst_recallprecision_plot[[dataset]], PCE == pce)
    p = ggplot(data, aes(y = precision, x = recall , color = method)) + 
      geom_point(size = 1.5) + 
      geom_line() + 
      geom_text_repel(aes(label = FC), size = 2.25,segment.linetype = 5,nudge_x = 0.005/max(data$recall)) + ggtitle(paste0("PCE = ",pce)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))  + 
      theme(legend.title = element_text(size = 7), 
            legend.text = element_text(size = 7),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),  
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            plot.title = element_text(size=10))
    lst_precision_recall_byPCE[[paste0("PCE_",pce)]] = p
  }

  ##################################
  ##### Generate f1 score plot #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  p1 = ggplot(tmp_df, aes(x = PCE, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~FC + PCE) +
    ggtitle("Faceted by FC")
  
  p2 = ggplot(tmp_df, aes(x = FC, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~PCE) +
    ggtitle("Faceted by PCE")
  
  ##########################################
  ##### Generate ranking LR genes plot #####
  
  tmp_df = ranking_LRgenes_lst_plot[[dataset]]
  p3 = ggplot(tmp_df, aes(x = FC, y = ratio, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~PCE) +
    ggtitle("Faceted by PCE . Ratio - n_simulated_LR in top50_significant_LR")
  
  ##### Save plots
  pdf(file.path(path_results_dir ,paste0(dataset, "_results_plots.pdf")), width = 12, height = 7)
  ggarrange(plotlist = lst_precision_recall_byPCE, common.legend = T) %>% print
  p1 %>% print
  p2 %>% print
  p3 %>% print
  heatmap1 %>% print
  heatmap2 %>% print
  
  ggarrange(plotlist = UpSet_plot_lst$VASAseq$FC_1_PCE_10,nrow = 2) %>% print
  ggarrange(plotlist = UpSet_plot_lst$VASAseq$FC_1_PCE_50,nrow = 2) %>% print
  ggarrange(plotlist = UpSet_plot_lst$VASAseq$FC_3_PCE_10,nrow = 2) %>% print
  ggarrange(plotlist = UpSet_plot_lst$VASAseq$FC_3_PCE_50,nrow = 2) %>% print
  ggarrange(plotlist = UpSet_plot_lst$VASAseq$FC_5_PCE_10,nrow = 2) %>% print
  ggarrange(plotlist = UpSet_plot_lst$VASAseq$FC_5_PCE_50,nrow = 2) %>% print
  
  dev.off()
}
