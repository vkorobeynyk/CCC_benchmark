# Load package
library(dplyr)
library(jsonlite)
library(stringr)

# An useful error if the argument is missing
if (is.null(snakemake@input[["significant_interactions"]]) | is.null(snakemake@input[["simulated_interactions"]]) | is.null(snakemake@output[["CT_statistics"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# OUTPUT FILES
path_CT_statistics <- snakemake@output[["CT_statistics"]]

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
        # create vectors with Interaction partners
        simulated_interactions_CTCT = simulated_interactions[[comb_CTs]]
        # remove the subunit genes that we inflated -> if we dont remove them, we will have much more FN
        simulated_interactions_CTCT = simulated_interactions_CTCT[!grepl("subunit", simulated_interactions_CTCT)]
        significant_interactions_CTCT = significant_interactions[[comb_CTs]] %>% mutate(LR = str_c(ligand, "_", receptor))

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

precision = lst_score_perCTCT[[1]]$precision%>% round(2)
recall =  lst_score_perCTCT[[1]]$recall %>% round(2)
f1score =  lst_score_perCTCT[[1]]$f1score %>% round(2)

# print all scores
print(paste("precision:", precision %>% round(4), "recall:" ,recall %>% round(4) , "f1score:", f1score %>% round(4)))

# Create df to store average results
df_statistics = data.frame(precision = precision, recall = recall, f1score = f1score)

# in case f1score is NaN
if(is.nan(df_statistics$f1score)) {df_statistics$f1score = 0}

######################
##### Save files #####
######################

write.csv(df_statistics, path_CT_statistics, row.names = F)