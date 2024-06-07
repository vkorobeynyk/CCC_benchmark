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

diagnostic_plots_realFC = list()
diagnostic_df_realFC = list()
diagnostic_plots_MeanVar = list()
for(dataset in datasets)
{
  file_path = file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB"))
  realFC_aftersimulation_files = file_path %>% list.files(., pattern = "realFC_aftersimulation")
  
  # Generate mean var plot of genes
  genemetadata = readRDS(file.path("data/processed/", dataset,"genemetadata.RDS"))
  p = ggplot(genemetadata$mean, aes(x = vars, y = means)) + geom_hex() +
    geom_point(data=genemetadata$mean[genemetadata$mean$gene_to_use,], aes(x=vars, y=means), colour="red", size=2) +
    xlab("variance (log10)") +
    ylab("mean (log10)") +
    ggtitle(paste0("Highlighted LR genes for ",dataset))
  
  diagnostic_plots_MeanVar[[dataset]] = p
  
  params_grid = expand.grid(vector1 = FC, vector2 = PCE)
  for(i in 1:nrow(params_grid))
  {
    x = params_grid[i,] %>% as.numeric ; names(x) = c("FC","PCE")
    
    ###########################################################
    ##### Generate Plots of real FC after semi-simulation #####
    ###########################################################
    
    realFC_aftersimulation_file = realFC_aftersimulation_files[grepl(paste0("^" , x["FC"] , "$"), str_split(realFC_aftersimulation_files, "_") %>% 
                                                                       lapply(., "[[", 3)) & grepl(paste0("^" , x["PCE"] , ".RDS$"), str_split(realFC_aftersimulation_files, "_") %>% lapply(., "[[", 5))]
    
    tmp_lst = readRDS(file.path(file_path,realFC_aftersimulation_file)) %>% 
      unlist %>% subset(.,!is.infinite(.)) %>% 
      plot_FCafter_semisimulation(. , theoreticalFC = x["FC"], PCE = x["PCE"])
    
    diagnostic_plots_realFC[[dataset]][[paste0("FC_",x["FC"])]][[paste0("PCE_",x["PCE"])]] = tmp_lst %>% pluck("plot")
    
    diagnostic_df_realFC[[dataset]][[paste0("index_",i)]] = tmp_lst[c(2,3,4)] %>% as.data.frame
  }
}

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
  #realFC_aftersimulation_files = file_path %>% list.files(., pattern = "realFC_aftersimulation")
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
    counts_file = inflated_counts_files[grepl(paste0("^" , x["FC"] , "$"), str_split(inflated_counts_files, "_") %>% lapply(., "[[", 5)) & 
                                          grepl(paste0("^" , x["PCE"] , ".tsv$"), str_split(inflated_counts_files, "_") %>% lapply(., "[[", 7))]
    master_lst_diagnosticPlots[[naming]][["counts"]] = read.table(file.path(file_path,counts_file))
    
    # load correct simulated interactions file depending on the params_grid
    simulated_interactions_file = simulated_interactions_files[grepl(paste0("^" , x["FC"] , "$"), str_split(simulated_interactions_files, "_") %>% 
                                                                       lapply(., "[[", 4)) & grepl(paste0("^" , x["PCE"] , ".RDS$"), str_split(simulated_interactions_files, "_") %>% lapply(., "[[", 6))]
    master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] = readRDS(file.path(file_path,simulated_interactions_file))
    
    ########################################
    ##### Generate L/R inflated per CT #####
    ########################################
    CT_present = names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_")  %>% unlist %>% unique
    for(CT in CT_present)
    {
      # Find Ligand genes which were inflated in specified celltype 
      n = which(names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_") %>% lapply(.,"[[",1) %>% unlist %in% CT)
      L_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]][n] %>% unlist %>% str_split(.,"_") %>% 
        lapply("[[", 1) %>% unlist %>% setdiff(., "subunit") # remove subunit string
      
      # Find Receptor genes which were inflated in specified celltype inflated 
      n = which(names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_") %>% lapply(.,"[[",2) %>% unlist %in% CT)
      R_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]][n] %>% unlist %>% str_split(.,"_") %>% 
        lapply("[[", 2) %>% unlist%>% setdiff(., "subunit") # remove subunit string
      
      # Save inflated genes per CT
      master_lst_diagnosticPlots[[naming]][[CT]] = list(L = L_genes_inflated_per_CT_to_keep %>% unique , R = R_genes_inflated_per_CT_to_keep %>% unique)
    }
    
    # load correct PCE file depending on the params_grid
    PCE_file = PCE_files[grepl(paste0("^" , x["FC"] , "$"), str_split(PCE_files, "_") %>% lapply(., "[[", 6)) & grepl(paste0("^" , x["PCE"] , ".RDS$"), str_split(PCE_files, "_") %>% lapply(., "[[", 8))]
    master_lst_diagnosticPlots[[naming]][["PCE"]] = readRDS(file.path(file_path,PCE_file)) %>% unlist * 100 # transform to percentage
  }
  # filter original avelogcpm counts to contain same genes as the inflated count matrices
  original_counts = original_counts %>% subset(rownames(original_counts) %in% (master_lst_diagnosticPlots[[1]]$counts %>% rownames))
  
  diagnostic_plots_lst[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
                                                             FC_param = FC, PCE_param = PCE, dataset = dataset, metadata = metadata, CT_toPlot = CT_present)
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
                                                               FC_param = FC, PCE_param = PCE, dataset = dataset, metadata = metadata, CT_toPlot = CT_present[1])
  #diagnostic_plots_lst[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
  #                                                           FC_param = FC[c(1,4,8)], PCE_param = PCE[c(1,3,7)], dataset = dataset, metadata = metadata, CT_toPlot = CT_present)
  #diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
  #                                                             FC_param = FC[c(1,4,8)], PCE_param = PCE[c(1,3,7)], dataset = dataset, metadata = metadata, CT_toPlot = CT_present[1])
}

############################### Diagnostic plots
for(dataset in datasets)
{
  pdf(file.path(path_results_dir ,paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  do.call(ggarrange,c(diagnostic_plots_lst[[dataset]]$avelogcpm_fixedPCE, common.legend = TRUE)) %>% print
  do.call(ggarrange,diagnostic_plots_lst[[dataset]]$PCE_fixedPCE) %>% print
  do.call(ggarrange,c(diagnostic_plots_perCT[[dataset]]$avelogcpm_fixedPCE, common.legend = TRUE)) %>% print
  do.call(ggarrange,diagnostic_plots_perCT[[dataset]]$PCE_fixedPCE) %>% print
  
  x = 1:length(diagnostic_plots_realFC[[dataset]])
  sapply(x, function(x) {do.call(ggarrange,diagnostic_plots_realFC[[dataset]][[x]]) %>% print}) %>% print
  diagnostic_plots_MeanVar[[dataset]] %>% print
  dev.off()
}

##################################
##### precision/recall plots #####
##################################

# plot TPR/sensitivity/recall
#statistics_results_lst = list()
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
      master_lst_precision_recall[[file]] = master_lst_precision_recall[[file]][c(1,2,3)] # remove the f1score
      
    }
    
    ###############################################################################
    ##### Ranking of LR genes across top 25% of significant hits from methods #####
    ###############################################################################
    tmp_ranking_LRgenes_lst = list()
    ranking_LRgenes_files = file.path(path_output_dir,dataset,"metrics/",method) %>%
      list.files(., pattern = "ranking_LRgenes")
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
                         select("FC_real_median") %>% as.numeric) %>% round
    }
    rm(tmp_diagnostic_df_realFC)
    
    statistics_results_lst_recallprecision_plot[[dataset]][[method]]$FC = FC_real
    statistics_results_lst_recallprecision_plot[[dataset]][[method]]$PCE = PCE
    statistics_results_lst_recallprecision_plot[[dataset]][[method]]$method = method
  }

  # for precision recall plots
  tmp = do.call(rbind, statistics_results_lst_recallprecision_plot[[dataset]])
  tmp$dataset = dataset
  statistics_results_lst_recallprecision_plot[[dataset]] = tmp
  
  # Add FC and PCE info for ranking LR genes plot 
  ranking_LRgenes_lst_plot[[dataset]] = ranking_LRgenes_lst_plot[[dataset]] %>% do.call(rbind,.) %>%
    mutate(PCE = statistics_results_lst_recallprecision_plot[[dataset]]$PCE ,
           FC = statistics_results_lst_recallprecision_plot[[dataset]]$FC)

  ##########################################
  ##### Generate Precision recall plot #####
  ##########################################
  
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
  ##################################
  
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
  ##########################################
  
  tmp_df = ranking_LRgenes_lst_plot[[dataset]]
  p3 = ggplot(tmp_df, aes(x = FC, y = ratio, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~PCE) +
    ggtitle("Faceted by PCE . Ration - n_simulated_LR in top25%_significant_LR")
  
  ##### Save plots
  pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots.pdf")), width = 12, height = 7)
  ggarrange(plotlist = lst_precision_recall_byPCE, common.legend = T) %>% print
  p1 %>% print
  p2 %>% print
  p3 %>% print
  dev.off()
}
