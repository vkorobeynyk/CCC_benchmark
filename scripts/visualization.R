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
library(cetcolor)
source("scripts/helper_functions.R")

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
path_config.yaml <- opt$config.yaml
path_output_dir <- opt$path_output_dir
path_results_dir <- opt$path_results_dir
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

#methods = methods[1:2]
#datasets = datasets[1:2]
#FC = FC[1:2]
#PCE = PCE[1:2]
############################
##### Diagnostic plots #####
############################

diagnostic_plots_lst = list()
diagnostic_plots_perCT= list()
for(dataset in datasets)
{
  
  # load files
  file_path = file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB"))
  inflated_counts_files = file_path %>% list.files(., pattern = "sc_inflated_counts")
  original_counts = read.table(file.path("data/",dataset, "raw_counts.tsv"))  # load original counts
  rownames(original_counts) = rownames(original_counts) %>% toupper()
  PCE_files = file_path %>% list.files(., pattern = "cells_sampled_perCTCT")
  metadata_files = file_path %>% list.files(., pattern = "metadata")
  simulated_interactions_files = file_path %>% list.files(., pattern = "simulated_interactions")
  
  # list to save all data
  master_lst_diagnosticPlots = list()
  
  # load metadata files
  # as the metadata files are the same -> load 1st one
  metadata_file = read.table(file.path("data/",dataset, "raw_metadata.tsv"))
  metadata_file$Celltype = gsub(" ", ".", metadata_file$Celltype)
  
  # generate the grid of parameters used for naming the list to generate outputs
  params_grid = expand.grid(vector1 = FC,
              vector2 = PCE)
  
  # iterate over the grid of parameters
  for(i in 1:nrow(params_grid))
  {
    x = params_grid[i,]
    
    naming = paste0("FC_",x[1],"_PCE_",x[2])
    
    # load counts
    counts_file = inflated_counts_files[grepl(paste0("^" , x[1] , "$"), str_split(inflated_counts_files, "_") %>% lapply(., "[[", 5)) & grepl(paste0("^" , x[2] , ".tsv$"), str_split(inflated_counts_files, "_") %>% lapply(., "[[", 7))]
    master_lst_diagnosticPlots[[naming]][["counts"]] = read.table(file.path(file_path,counts_file))
    gene_names = rownames(master_lst_diagnosticPlots[[naming]][["counts"]])
    
    # load simulated interactions
    simulated_interactions_file = simulated_interactions_files[grepl(paste0("^" , x[1] , "$"), str_split(simulated_interactions_files, "_") %>% lapply(., "[[", 4)) & grepl(paste0("^" , x[2] , ".RDS$"), str_split(simulated_interactions_files, "_") %>% lapply(., "[[", 6))]
    master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] = readRDS(file.path(file_path,simulated_interactions_file))
    
    ########################################
    ##### Generate L/R inflated per CT #####
    ########################################
    # metadata originally is not 
    CT_present = which(unique(metadata_file$Celltype) %in% (names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_") %>% lapply(.,"[[",1) %>% unlist %>% unique))
    for(CT in unique(metadata_file$Celltype)[CT_present])
    {
      # Find Ligand genes which were inflated in specified celltype inflated 
      n = which(names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_") %>% lapply(.,"[[",1) %>% unlist %in% CT)
      L_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]][n] %>% unlist %>% str_split(.,"_") %>% lapply("[[", 1) %>% unlist
      L_genes_inflated_per_CT_to_keep = L_genes_inflated_per_CT_to_keep[!grepl("subunit", L_genes_inflated_per_CT_to_keep)] # remove subunit string
      
      # Find Receptor genes which were inflated in specified celltype inflated 
      n = which(names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_") %>% lapply(.,"[[",2) %>% unlist %in% CT)
      R_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]][n] %>% unlist %>% str_split(.,"_") %>% lapply("[[", 2) %>% unlist
      R_genes_inflated_per_CT_to_keep = R_genes_inflated_per_CT_to_keep[!grepl("subunit", R_genes_inflated_per_CT_to_keep)] # remove subunit string
      
      # Save inflated genes per CT
      master_lst_diagnosticPlots[[naming]][[CT]] = list(L = L_genes_inflated_per_CT_to_keep %>% unique , R = R_genes_inflated_per_CT_to_keep %>% unique)
    }
    
    # load PCE information
    PCE_file = PCE_files[grepl(paste0("^" , x[1] , "$"), str_split(PCE_files, "_") %>% lapply(., "[[", 5)) & grepl(paste0("^" , x[2] , ".RDS$"), str_split(PCE_files, "_") %>% lapply(., "[[", 7))]
    master_lst_diagnosticPlots[[naming]][["PCE"]] = readRDS(file.path(file_path,PCE_file)) %>% unlist * 100 # transform to percentage
  }
  # fitler original avelogcpm counts to contain same genes as the inflated count matrices
  original_counts = original_counts %>% subset(rownames(original_counts) %in% (master_lst_diagnosticPlots[[1]]$counts %>% rownames))
  
  #diagnostic_plots_lst[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
  #                                                           FC_param = FC, PCE_param = PCE, dataset = dataset, cell_metadata = metadata_file , CT_toPlot = unique(metadata_file$Celltype))
  #diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
  #                                                           FC_param = FC, PCE_param = PCE, dataset = dataset, cell_metadata = metadata_file , CT_toPlot = unique(metadata_file$Celltype)[1])
  diagnostic_plots_lst[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
                                                             FC_param = FC[c(1,4,8)], PCE_param = PCE[c(1,4,8)], dataset = dataset, cell_metadata = metadata_file , CT_toPlot = unique(metadata_file$Celltype))
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
                                                               FC_param = FC[c(1,4,8)], PCE_param = PCE[c(1,4,8)], dataset = dataset, cell_metadata = metadata_file , CT_toPlot = unique(metadata_file$Celltype)[1])
  
}


############################### Diagnostic plots
for(dataset in datasets)
{
  pdf(file.path(path_results_dir ,paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  do.call(ggarrange,diagnostic_plots_lst[[dataset]]$avelogcpm_fixedPCE) %>% print
  do.call(ggarrange,diagnostic_plots_lst[[dataset]]$PCE_fixedPCE) %>% print
  do.call(ggarrange,diagnostic_plots_perCT[[dataset]]$avelogcpm_fixedPCE) %>% print
  do.call(ggarrange,diagnostic_plots_perCT[[dataset]]$PCE_fixedPCE) %>% print
  dev.off()
}

##################################
##### precision/recall plots #####
##################################

# plot TPR/sensitivity/recall
averaged_statistics_results_lst = list()
averaged_statistics_results_lst_recallprecision_plot = list()
master_lst_precision_recall = list()
for(dataset in datasets)
{
  for(method in methods)
  {
    averaged_statistics_files = file.path(path_output_dir,dataset,"metrics/",method) %>%
      list.files(., pattern = "averaged_statistics")
    master_lst_precision_recall = list()
    for(file in averaged_statistics_files)
    {
      master_lst_precision_recall[[file]] = read.csv((file.path(path_output_dir,dataset,"metrics/",method,file))) %>% unlist
      master_lst_precision_recall[[file]] = master_lst_precision_recall[[file]][c(1,2,4,5)] # remove the f1score
      
    }
    # Generate data for plots isolation FC/PCE 
    x1 = melt(do.call(rbind,master_lst_precision_recall)[,c(1,2)])
    colnames(x1) =  c("name", "metric" , "value")
    x2 = melt(do.call(rbind,master_lst_precision_recall)[,c(3,4)])
    colnames(x2) =  c("name", "metric_sd" , "value_sd")
    averaged_statistics_results_lst[[dataset]][[method]] = cbind(x1,x2)
    averaged_statistics_results_lst[[dataset]][[method]] = averaged_statistics_results_lst[[dataset]][[method]][,c(1,2,3,5,6)] # remove extra name column
    
    # Add FC and PCE columns
    x = averaged_statistics_results_lst[[dataset]][[method]]$name %>% stringr::str_split("_")
    tmp_FC = lapply(x, "[[" , 5) %>% unlist %>% as.double
    PCE = lapply(x, "[[" , 7) %>% unlist %>% as.double
    
    averaged_statistics_results_lst[[dataset]][[method]]$FC = tmp_FC
    averaged_statistics_results_lst[[dataset]][[method]]$PCE = PCE
    averaged_statistics_results_lst[[dataset]][[method]]$method = method
    
    # Generate data for precision vs recall plots
    averaged_statistics_results_lst_recallprecision_plot[[dataset]][[method]] = do.call(rbind,master_lst_precision_recall) %>% as.data.frame
    
    # Add FC and PCE columns
    x = rownames(averaged_statistics_results_lst_recallprecision_plot[[dataset]][[method]]) %>% stringr::str_split("_")
    tmp_FC = lapply(x, "[[" , 5) %>% unlist %>% as.double
    PCE = lapply(x, "[[" , 7) %>% unlist %>% as.double
    
    averaged_statistics_results_lst_recallprecision_plot[[dataset]][[method]]$FC = tmp_FC
    averaged_statistics_results_lst_recallprecision_plot[[dataset]][[method]]$PCE = PCE
    averaged_statistics_results_lst_recallprecision_plot[[dataset]][[method]]$method = method
  }
  
  tmp = do.call(rbind, averaged_statistics_results_lst[[dataset]])
  tmp$dataset = dataset
  averaged_statistics_results_lst[[dataset]] = tmp
  
  tmp = do.call(rbind, averaged_statistics_results_lst_recallprecision_plot[[dataset]])
  tmp$dataset = dataset
  averaged_statistics_results_lst_recallprecision_plot[[dataset]] = tmp
  
  
  ##########################################
  ##### Generate Precision recall plot #####
  ##########################################
  
  PCE = as.factor(averaged_statistics_results_lst_recallprecision_plot[[dataset]]$PCE)
  FC = as.factor(averaged_statistics_results_lst_recallprecision_plot[[dataset]]$FC)
  
  # PLot precision recall curves when generating 1 plot per FC
  lst_precision_recall_byFC = list()
  for(fc in levels(FC))
  {
    data = filter(averaged_statistics_results_lst_recallprecision_plot[[dataset]], FC == fc)
    p = ggplot(data, aes(y = averaged_precision, x = averaged_recall , color = method)) + 
      geom_point(size = 1.5) + 
      geom_line() + 
      geom_text_repel(aes(label = PCE), size = 2.25) + ggtitle(paste0("FC = ",fc)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
      theme(legend.title = element_text(size = 7), 
            legend.text = element_text(size = 7),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),  
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            plot.title = element_text(size=10))
    lst_precision_recall_byFC[[paste0("FC_",fc)]] = p
  }
  
  # PLot precision recall curves when generating 1 plot per PCE
  lst_precision_recall_byPCE = list()
  for(pce in levels(PCE))
  {
    data = filter(averaged_statistics_results_lst_recallprecision_plot[[dataset]], PCE == pce)
    p = ggplot(data, aes(y = averaged_precision, x = averaged_recall , color = method)) + 
      geom_point(size = 1.5) + 
      geom_line() + 
      geom_text_repel(aes(label = FC), size = 2.25) + ggtitle(paste0("PCE = ",pce)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
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
  # Fixed FC
  lst_valueVSpce = list()
  for(fc in unique(FC) %>% sort)
  {
    data = filter(averaged_statistics_results_lst[[dataset]], FC == fc)
    p1 = ggplot(filter(data, metric == "averaged_precision"), aes(y = value, x = PCE , color = method)) + 
      geom_line() + 
      ggtitle(paste("FC =" , fc, "dataset" , dataset)) + facet_grid(~metric) #  + 
      # geom_label_repel(aes(label = method) ,nudge_x= 0.25, segment.size= 0.2, max.overlaps = 30) + guides(color= guide_legend(override.aes = aes(label = "")))
    
    p2 = ggplot(filter(data, metric == "averaged_recall"), aes(y = value, x = PCE , color = method)) + 
      geom_line() + 
      ggtitle(paste("FC =" , fc, "dataset" , dataset)) + facet_grid(~metric) #+ geom_label_repel(aes(label = method) ,nudge_x= 0.25, segment.size= 0.2)
    
    lst_valueVSpce[[paste0("FC_",fc)]] = ggarrange(p1,p2, common.legend = T)
  }
  
  # Fixed % cells expressing
  lst_valueVSfc = list()
  for(pce in unique(PCE) %>% sort)
  {
    data = filter(averaged_statistics_results_lst[[dataset]], PCE == pce)
    p1 = ggplot(filter(data, metric == "averaged_precision"), aes(y = value, x = FC , color = method)) + 
      geom_line() + 
      ggtitle(paste("PCE =" , pce, "dataset" , dataset)) + facet_grid(~metric) #+ geom_label_repel(aes(label = method) ,nudge_x= 0.25, segment.size= 0.2)
    p2 = ggplot(filter(data, metric == "averaged_recall"), aes(y = value, x = FC , color = method)) + 
      geom_line() + 
      ggtitle(paste("PCE =" , pce, "dataset" , dataset)) + facet_grid(~metric) #+ geom_label_repel(aes(label = method) ,nudge_x= 0.25, segment.size= 0.2)
    
    lst_valueVSfc[[paste0("PCE",pce)]] = ggarrange(p1,p2, common.legend = T)
  }
  
  
  pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots.pdf")), width = 12, height = 7)
  ggarrange(plotlist = lst_precision_recall_byFC, common.legend = T) %>% print
  ggarrange(plotlist = lst_precision_recall_byPCE, common.legend = T) %>% print
  lst_valueVSpce %>% print
  lst_valueVSfc %>% print
  dev.off()
}
