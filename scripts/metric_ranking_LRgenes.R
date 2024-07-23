# Load package
library(dplyr)
library(jsonlite)
library(stringr)

# An useful error if the argument is missing
if (is.null(snakemake@input[["significant_interactions"]]) | is.null(snakemake@input[["simulated_interactions"]]) | is.null(snakemake@output[["ranking_LRgenes"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# OUTPUT FILES
path_ranking_LRgenes <- snakemake@output[["ranking_LRgenes"]]

# INPUT FILES
path_significant_interactions <- snakemake@input[["significant_interactions"]]
path_simulated_interactions <- snakemake@input[["simulated_interactions"]]

########################################
##### Loading and processing files #####
########################################

# load significant interactions
significant_interactions = readRDS(path_significant_interactions)

# load simulated interactions (TRUE)
simulated_interactions = readRDS(path_simulated_interactions)


# replace weird symbols to underscore in list names
names(significant_interactions) = significant_interactions %>% names %>% gsub( "[|]" , "_" , .)


#############################################
##### Generating a ranking for LR genes #####
#############################################
# By ranking I mean select 25% of genes with highest logFC (for example) and then check how many of those pairs we simulated

lst_score_perCTCT = list()
# extract all possible CTs combinations
combinations_CTs = names(simulated_interactions)

# loop over all CT combination
for(comb_CTs in combinations_CTs) 
{
  if(comb_CTs %in% names(significant_interactions))
  {
    # some method may have no significant interaction at lower parameter values -> skip those
    if(dim(significant_interactions[[comb_CTs]])[1] == 0) {
      df_statistics = data.frame(amount_simulatedLR_intop25perc = 0,
                                 amount_significantLR_intop25perc = 0,
                                 ratio = 0)
      next
    }
    # create vectors with Interaction partners
    simulated_interactions_CTCT = simulated_interactions[[comb_CTs]]
    # remove the subunit genes that we inflated -> if we dont remove them, we will have much more FN
    simulated_interactions_CTCT = simulated_interactions_CTCT[!grepl("subunit", simulated_interactions_CTCT)]
    
    #####################################################################################################################
    ### Almost every method has different ways to compute strength of interaction -> so this part is always different ###
    #####################################################################################################################
    if(significant_interactions[[comb_CTs]]$method[1] == "cellphonedb")
    {
      significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR_logFC = abs(LR_logFC) , LR = str_c(ligand, "_", receptor)) %>% 
        arrange(.,desc(LR_logFC)) 
      
      # select 25% genes with highest logFC for ranking
      significant_interactions_CTCT = significant_interactions_CTCT %>% top_n(nrow(significant_interactions_CTCT) * 0.25)
    } else if(significant_interactions[[comb_CTs]]$method[1] == "cellchat")
    {
      significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR = str_c(ligand, "_", receptor)) %>% 
        arrange(.,desc(prob)) 
      
      significant_interactions_CTCT = significant_interactions_CTCT %>% top_n(nrow(significant_interactions_CTCT) * 0.25)
    } else if(significant_interactions[[comb_CTs]]$method[1] == "connectome")
    {
      significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR = str_c(ligand, "_", receptor)) %>% 
        arrange(.,desc(weight_sc)) # we are adding expression so we should see increase in specificity == increase in weight_sc
      
      significant_interactions_CTCT = significant_interactions_CTCT %>% top_n(nrow(significant_interactions_CTCT) * 0.25)
    } else if(significant_interactions[[comb_CTs]]$method[1] == "cytotalk")
    {
      significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR = str_c(ligand, "_", receptor)) %>% 
        arrange(.,desc(crosstalk_score)) 
      
      significant_interactions_CTCT = significant_interactions_CTCT %>% top_n(nrow(significant_interactions_CTCT) * 0.25)
    } else if(significant_interactions[[comb_CTs]]$method[1] == "singlecellsignalR")
    {
      significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR = str_c(ligand, "_", receptor)) %>% 
        arrange(.,desc(LRscore)) 
      
      # select 25% genes with highest logFC for ranking
      significant_interactions_CTCT = significant_interactions_CTCT %>% top_n(nrow(significant_interactions_CTCT) * 0.25)
    } else if(significant_interactions[[comb_CTs]]$method[1] == "natmi_specificity")
    {
      significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR = str_c(ligand, "_", receptor)) %>% 
        arrange(.,desc(edge_specificity)) 
      
      # select 25% genes with highest logFC for ranking
      significant_interactions_CTCT = significant_interactions_CTCT %>% top_n(nrow(significant_interactions_CTCT) * 0.25)
    } else if(significant_interactions[[comb_CTs]]$method[1] == "natmi_4th_quantile")
    {
      significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR = str_c(ligand, "_", receptor)) %>% 
        arrange(.,desc(sum_LR_specificity)) 
      
      # select 25% genes with highest logFC for ranking
      significant_interactions_CTCT = significant_interactions_CTCT %>% top_n(nrow(significant_interactions_CTCT) * 0.25)
    } else{message(paste("the method",significant_interactions[[comb_CTs]]$method[1], "was not accounted for in ranking_LRgenes script" ))}
    
    
    
    # how many simulated gene pairs are in the methods output
    df_statistics = data.frame(amount_simulatedLR_intop25perc = intersect(simulated_interactions_CTCT, significant_interactions_CTCT$LR) %>% length,
                    amount_significantLR_intop25perc = nrow(significant_interactions_CTCT)) %>% mutate(ratio = amount_simulatedLR_intop25perc / amount_significantLR_intop25perc)
  } 
}

######################
##### Save files #####
######################

write.csv(df_statistics, path_ranking_LRgenes, row.names = F)