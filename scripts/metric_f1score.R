# Load package
library(dplyr)
library(jsonlite)

# An useful error if the argument is missing
if (is.null(snakemake@input[["significant_interactions"]]) | is.null(snakemake@input[["simulated_interactions"]]) | is.null(snakemake@output[["averaged_statistics"]]) | is.null(snakemake@output[["CT_statistics"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read the argument
path_significant_interactions <- snakemake@input[["significant_interactions"]]
path_simulated_interactions <- snakemake@input[["simulated_interactions"]]
path_averaged_statistics <- snakemake@output[["averaged_statistics"]]
path_CT_statistics <- snakemake@output[["CT_statistics"]]


########################################
##### Loading and processing files #####
########################################

# load significant interactions
significant_interactions = readRDS(path_significant_interactions)

# load simulated interactions (TRUE)
simulated_interactions = readRDS(path_simulated_interactions)


# replace weird symbols to underscore in list names
names(significant_interactions) = significant_interactions %>% names %>% gsub( "[|]" , "_" , .)


####################################################
##### Generating precision,recall and f1scores #####
####################################################

lst_score_perCTCT = list()
# extract all possible CTs combinations
combinations_CTs = names(simulated_interactions)

# loop over all CT combination
for(comb_CTs in combinations_CTs) 
{
    if(comb_CTs %in% names(significant_interactions))
    {
        # creave vectors with Interaction partners
        simulated_interactions_CTCT = simulated_interactions[[comb_CTs]]
        # remove the subunit genes that we inflated -> if we dont remove them, we will have much more FN
        simulated_interactions_CTCT = simulated_interactions_CTCT[!grepl("subunit", simulated_interactions_CTCT)]
        significant_interactions_CTCT = significant_interactions[[comb_CTs]]
        significant_interactions_CTCT$LR = paste0(significant_interactions_CTCT$ligand , "_" , significant_interactions_CTCT$receptor)

        # compute statistics
        TP = intersect(simulated_interactions_CTCT, significant_interactions_CTCT$LR) %>% length
        #TN = 
        FP = setdiff(significant_interactions_CTCT$LR , simulated_interactions_CTCT) %>% length
        FN = setdiff(simulated_interactions_CTCT , significant_interactions_CTCT$LR) %>% length

        if(TP == 0) {precision = 0 ; recall = 0} else {precision = TP / (TP + FP) ; recall = TP / (TP + FN)}
        
        
        lst_score_perCTCT[[comb_CTs]][["precision"]] = precision %>% round(2)
        lst_score_perCTCT[[comb_CTs]][["recall"]] = recall %>% round(2)
        lst_score_perCTCT[[comb_CTs]][["f1score"]] = 2 * precision * recall / (precision + recall) %>% round(2)
    } else {
        lst_score_perCTCT[[comb_CTs]][["precision"]] = NA
        lst_score_perCTCT[[comb_CTs]][["recall"]] = NA
        lst_score_perCTCT[[comb_CTs]][["f1score"]] = NA
    }
    
    # In case precision and recall are 0 , replace NaN by 0
    if(is.nan(lst_score_perCTCT[[comb_CTs]][["f1score"]])) {lst_score_perCTCT[[comb_CTs]][["f1score"]] = 0}
}

##############
# statistics #
##############

# calculate average statistics across all CT_CT combination
# MEAN
mean_precision = sapply(combinations_CTs, function(x) { lst_score_perCTCT[[x]]$precision }) %>% mean %>% round(2)
mean_recall = sapply(combinations_CTs, function(x) { lst_score_perCTCT[[x]]$recall }) %>% mean %>% round(2)
mean_f1score = sapply(combinations_CTs, function(x) 
{ 
    # sometimes precision and recall is 0 which gives NaN f1score -> replace by 0
    if(lst_score_perCTCT[[x]]$f1score %>% is.nan) {0} else {lst_score_perCTCT[[x]]$f1score }#

}) %>% mean %>% round(2)

# create df to store results per CT_CT combination
df_score_perCTcombinations = do.call(rbind.data.frame, lst_score_perCTCT)
df_score_perCTcombinations$interacting_CT = rownames(df_score_perCTcombinations)

# SD
sd_precision = sd(df_score_perCTcombinations$precision) %>% round(3)
sd_recall = sd(df_score_perCTcombinations$recall) %>% round(3)
sd_f1score = sd(df_score_perCTcombinations$f1score) %>% round(3)

# print all scores
print(paste("averaged precision:", mean_precision %>% round(4), "averaged recall:" ,mean_recall %>% round(4) , "averaged f1score:", mean_f1score %>% round(4)))

# Create df to store average results
df_statistics = data.frame(averaged_precision = mean_precision, averaged_recall = mean_recall, averaged_f1score = mean_f1score,
                           sd_precision = sd_precision, sd_recall = sd_recall, sd_f1score = sd_f1score)

# in case f1score is NaN
if(is.nan(df_statistics$averaged_f1score)) {df_statistics$averaged_f1score = 0}

######################
##### Save files #####
######################

write.csv(df_statistics, path_averaged_statistics, row.names = F)
write.csv(df_score_perCTcombinations, path_CT_statistics, row.names = F)
