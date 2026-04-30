suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(glue)
})

# An useful error if the argument is missing
if (is.null(snakemake@input[["significant_interactions"]]) |  is.null(snakemake@output[["final_scores"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# OUTPUT FILES
path_final_scores = snakemake@output[["final_scores"]]

# INPUT FILES
path_significant_interactions = snakemake@input[["significant_interactions"]] # not used in the script but only in Snakefile

##########
# PARAMS #
##########

config = yaml::read_yaml("config.yaml")

# Setting parameters
methods = config$methods %>% unlist
datasets = config$datasets %>% unlist
FC = config$semiSimulation$FC %>% unlist %>% as.double
PCE = config$semiSimulation$PCE %>% unlist %>% as.double

# df: the ordered output from one method
calculate_normalizedRank = function(df, target_genes) {
  
  N = nrow(df)
  # get ranks
  ranks = match(target_genes, df$ligand_receptor)
  # Handle not found genes
  ranks[is.na(ranks)] = N 
  # Normalize the ranks
  norm_ranks = (ranks - 1) / (N - 1)
  
  return(mean(norm_ranks))
}

########################################
##### Loading and processing files #####
########################################
statistics_lst = list()
for(dataset in datasets)
{
  # load simulated interactions -> same indexLR has same interaction
  simulated_interaction_file = list.files(paste0("output/",dataset,"_semiSimulation_NB"), pattern = "simulated_interactions")[1]
  simulated_interactions = readRDS(paste0("output/",dataset,"_semiSimulation_NB/",simulated_interaction_file ))
  
  for(method in methods)
  {
    # load all output files for every combination of parameter for specific method
    data_dir =  file.path("output/",dataset,method)
    
    # iterate over every combination of parameters and retrieve statistics
    params_grid = expand.grid(FC = FC, PCE = PCE) %>%
      mutate(file_path = glue("{data_dir}/significant_interactions_FC_{FC}_PercCellsExpressing_{PCE}.csv"))
    
    # iterate over the grid of parameters
    for(file in params_grid$file_path)
    {
      x = params_grid %>% filter(file_path == file)
      
      naming = paste0("FC_",x["FC"], "_PCE_", x["PCE"])
      
      # load correct count file depending on the params_grid
      all_interactions = read.table(file, header = T) 
      
      # select significant interactions}
      significant_interactions = all_interactions %>% filter(as.logical(significant))
      
      #############################################
      ##### Generating a ranking for LR genes #####
      #############################################
      
      perc_significantInteractions_among_simulated = mean(simulated_interactions$ligand_receptor %in% significant_interactions$ligand_receptor) * 100
      
      #  Get the ranks of genes from simulated_interactions found in significant_interactions
      normalizedRank_simulated_among_significant = calculate_normalizedRank(significant_interactions , target_genes = simulated_interactions$ligand_receptor)
      if(is.nan(normalizedRank_simulated_among_significant)) {normalizedRank_simulated_among_significant = 0}
      
      ####################################################
      ##### Generating precision,recall and f1scores #####
      ####################################################
      
      lst_score_perCTCT = list()
      
      TP = sum(simulated_interactions$ligand_receptor %in% significant_interactions$ligand_receptor) # LR that are in significant and simulated
      
      FP = sum(!significant_interactions$ligand_receptor %in% simulated_interactions$ligand_receptor) # LR that are in significant but not simulated
      
      FN = sum(!simulated_interactions$ligand_receptor %in% significant_interactions$ligand_receptor) # LR from simulated that are not in the significant
      
      if(TP == 0) {precision = 0 ; recall = 0} else {precision = TP / (TP + FP) ; recall = TP / (TP + FN)}
      f1score = 2 * precision * recall / (precision + recall) %>% round(3)
      if(is.nan(f1score)) {f1score = 0}
      
      # Create df to store average results
      df_statistics = data.frame(precision = precision %>% round(3), recall = recall %>% round(3), f1score = f1score %>% round(3), 
                                 TP = TP, FN = FN, FP = FP,
                                 method = method, 
                                 FC = x["FC"], 
                                 PCE = x["PCE"],
                                 perc_significantInteractions_among_simulated = perc_significantInteractions_among_simulated,
                                 normalizedRank_simulated_among_significant = normalizedRank_simulated_among_significant)
      
      statistics_lst[[dataset]][[method]][[naming]] = df_statistics
    }
  }
}

######################
##### Save files #####
######################

saveRDS(statistics_lst, path_final_scores)