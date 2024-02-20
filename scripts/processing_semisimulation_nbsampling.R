library(dplyr)
library(stringr)
library(liana)
library(magrittr)

# An useful error if the argument is missing
if (is.null(snakemake@input[["counts_processed"]]) | is.null(snakemake@input[["metadata_processed"]]) | is.null(snakemake@input[["gene_metadata"]]) | is.null(snakemake@input[["target_ct_file"]]) | 
    is.null(snakemake@params[["nLR_per_CTCTcomb"]]) | is.null(snakemake@wildcards[["FC"]]) | is.null(snakemake@wildcards[["perc_cells_expressing"]]) | 
    is.null(snakemake@output[["sc_inflated_counts"]]) | is.null(snakemake@output[["simulated_interactions"]])  | is.null(snakemake@output[["sc_metadata"]])  | is.null(snakemake@output[["cells_sampled_perCTCT"]])){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read the argument
path_sc_inflated_counts <- snakemake@input[["sc_inflated_counts"]]
path_sc_metadata <- snakemake@input[["sc_metadata"]]
path_significant_interactions <- snakemake@output[["significant_interactions"]]
path_simulated_interactions <- snakemake@input[["simulated_interactions"]]


# Call the argument
counts_processed_path <- snakemake@input[["counts_processed"]]
metadata_processed_path <- snakemake@input[["metadata_processed"]]
gene_metadata_path <- snakemake@input[["gene_metadata"]]
target_ct_file_path <- snakemake@input[["target_ct_file"]]
nLR_per_CTCTcomb <- snakemake@params[["nLR_per_CTCTcomb"]]
FC <- as.double(snakemake@wildcards[["FC"]])
perc_cells_expressing <- as.integer(snakemake@wildcards[["perc_cells_expressing"]])


# output files
sc_inflated_counts_path <- snakemake@output[["sc_inflated_counts"]]
simulated_interactions_path <- snakemake@output[["simulated_interactions"]]
sc_metadata_path <- snakemake@output[["sc_metadata"]]
cells_sampled_perCTCT_path <- snakemake@output[["cells_sampled_perCTCT"]]

##################
# Set parameters #
##################

nLR_per_CTCTcomb = nLR_per_CTCTcomb # LR pairs for each CT/CT combination
FC = FC # extent of count inflation (FC 2 = 2* mu)

#############
# read data #
#############

counts = read.table(counts_processed_path)
rownames(counts) = toupper(rownames(counts))
metadata = read.table(metadata_processed_path)
rownames(metadata) = metadata$cell_ID
genemetadata = readRDS(gene_metadata_path)
rownames(genemetadata$disp) = genemetadata$disp$gene %>% toupper
means_perCT = genemetadata$mean
rownames(means_perCT) = rownames(means_perCT) %>% toupper
target_ct = read.table(target_ct_file_path) %>% unlist %>% as.character

#counts = read.table("/home/vkorob/Documents/snakemake_CCC/benchmark_CCC_sc/data/processed/VASAseq/counts_VASAseq_processed.tsv")
#rownames(counts) = toupper(rownames(counts))
#metadata = read.table("/home/vkorob/Documents/snakemake_CCC/benchmark_CCC_sc/data/processed/VASAseq/metadata_VASAseq_processed.tsv")
#rownames(metadata) = metadata$cell_ID
#genemetadata = readRDS("/home/vkorob/Documents/snakemake_CCC/benchmark_CCC_sc/data/processed/VASAseq/gene_metadata.tsv")
#rownames(genemetadata$disp) = genemetadata$disp$gene %>% toupper
#means_perCT = genemetadata$mean
#rownames(means_perCT) = rownames(means_perCT) %>% toupper
#target_ct = read.table("/home/vkorob/Documents/snakemake_CCC/benchmark_CCC_sc/data/processed/VASAseq//target_ct_file.tsv") %>% unlist %>% as.character


# check if cell names of counts and metadata correspond and are in the same order
stopifnot(colnames(counts) == metadata$cell_ID)

# Subset counts and metadata to contain only the selected CT
n = which(metadata$Celltype %in% target_ct)
counts = counts[,n]
metadata = metadata[n,]

# remove genes
message(paste("Target celltypes:" , str_flatten(target_ct, " ")))
message(paste("Shape of count dataframe:" , str_flatten(dim(counts) , " ")))


#############
# Load LRdb #
#############

# Load OmniPath database
LRdb = select_resource(c('OmniPath'))[[1]]
colnames(LRdb)[1:2] = c("ligand", "receptor")

LRdb = LRdb[(LRdb$ligand %in% rownames(counts)),]
LRdb = LRdb[(LRdb$receptor %in% rownames(counts)),]
# filter LR because there are duplicated pairs (only 1)
LRdb$L_R = str_c(LRdb$ligand,"_",LRdb$receptor)
LRdb = LRdb[!duplicated(LRdb$L_R),]

message(paste("After filtering LR database," , nrow(LRdb) , "LR pairs show expression in at least 10 cells"))

###########################
# Inflate gene expression #
###########################
set.seed(3)
combination_CT = expand.grid(target_ct,target_ct)
combination_CT = paste0(combination_CT$Var1, "_", combination_CT$Var2)

simulated_interactions_lst = list()

tmp_LRdb = LRdb

# Pre-sample LR pairs to be used in the semi simulation
# Also pre-sample the subunit genes
for(comb_CT in combination_CT)
{
  # select randomly LR pair to inflate expression
  index_toSample = sample(seq(1,nrow(tmp_LRdb)) , size = nLR_per_CTCTcomb, replace = F)
  LR_sample = tmp_LRdb[index_toSample,]
  
  # Add genes belonging to a complex to be also inflated (based on CellPhoneDB and CellChatDB DB)
  # CellChat
  LRdb_cellchat = select_resource(c('CellChatDB'))[[1]]
  LRdb_cellchat = rbind(filter(LRdb_cellchat, grepl("COMPLEX", source)) , filter(LRdb_cellchat, grepl("COMPLEX", target)))
  
  # CellPhoneDB
  LRdb_cpdb = select_resource(c('CellPhoneDB'))[[1]]
  LRdb_cpdb = rbind(filter(LRdb_cpdb, grepl("COMPLEX", source)) , filter(LRdb_cpdb, grepl("COMPLEX", target)))
  
  LRdb_join = rbind(LRdb_cellchat , LRdb_cpdb)
  LRdb_join = LRdb_join[!str_c(LRdb_join$source, LRdb_join$target) %>% duplicated,] # remove duplicated entries
  
  # remove the sampled LR pairs 
  tmp_LRdb = tmp_LRdb[-index_toSample,]
  
  # adding L and R to the LR_list to track inflated genes without removing duplicates
  simulated_interactions_lst[[comb_CT]] = LR_sample$L_R
  
  ############## find subunits for ligands
  for(gene in LR_sample$ligand)
  {
    tmp_df = filter(LRdb_join , grepl(gene, source_genesymbol) & grepl("COMPLEX", source))
    
    # in case this gene has no subunits, skip
    if(nrow(tmp_df) == 0) {next}

    subunits = c(tmp_df$source_genesymbol %>% str_split("_") %>% lapply("[",1) , tmp_df$source_genesymbol %>% str_split("_") %>% lapply("[",2)) %>% unlist %>% unique()
    subunits = subunits[!grepl(gene, subunits)] # remove original gene
    
    simulated_interactions_lst[[comb_CT]]  %<>% append(. , subunits[which(subunits %in% rownames(counts))] %>% str_c(., "_subunit")) # remove empty strings and add subunit . Also here we filter subunits that are not present in count data
  }
  
  ############## find subunits for receptors
  for(gene in LR_sample$receptor)
  {
    tmp_df = filter(LRdb_join , grepl(gene, target_genesymbol) & grepl("COMPLEX", target))
    
    # in case this gene has no subunits, skip
    if(nrow(tmp_df) == 0) {next}
    
    subunits = c(tmp_df$target_genesymbol %>% str_split("_") %>% lapply("[",1) , tmp_df$target_genesymbol %>% str_split("_") %>% lapply("[",2)) %>% unlist %>% unique()
    subunits = subunits[!grepl(gene, subunits)] # remove original gene
    
    simulated_interactions_lst[[comb_CT]]  %<>% append(. , subunits[which(subunits %in% rownames(counts))] %>% str_c("subunit_" , .)) # remove empty strings and add subunit . Also here we filter subunits that are not present in count data
  }
  
  # remove duplicates in subunits info
  simulated_interactions_lst[[comb_CT]] = simulated_interactions_lst[[comb_CT]][!duplicated(simulated_interactions_lst[[comb_CT]])]
}

# Inflate expression of pre-sampled genes and in the respective CTs
inflated_LR_counts_lst = list()
perc_cells_expressing_lst = list()
counts_inflated = counts
rownames(LRdb) = LRdb$L_R

for(comb_CT in combination_CT)
{
  tmp_var1 = str_split(comb_CT,"_")[[1]]
  CT_sender = tmp_var1[1]
  CT_receiver = tmp_var1[2]

  L_sample = simulated_interactions_lst[[
    
  ]] %>% str_split("_") %>% lapply(.,"[[",1) %>% as.character
  R_sample = simulated_interactions_lst[[comb_CT]] %>% str_split("_") %>% lapply(.,"[[",2) %>% as.character
  
  # remove subunit string from the L and R vectors
  L_sample = L_sample[which(!L_sample %in% "subunit")]
  R_sample = R_sample[which(!R_sample %in% "subunit")]

    # iterate over cell type combination
  for(tmp_CT in c("CTsender","CTreceiver"))  
  {
    if (tmp_CT == "CTsender" ) {genes_to_sample = L_sample ; CT = CT_sender
    } else if (tmp_CT == "CTreceiver") {genes_to_sample = R_sample ; CT = CT_receiver}
    
    # select cells belonging to CT
    CT_cells = colnames(counts)[which(metadata$Celltype == CT)]
    # only select specific percentage of cells to increase expression
    cells_to_impute = sample(CT_cells, (length(CT_cells) * perc_cells_expressing / 100) %>% ceiling)
    # Iterate over every gene (L/R) depending on the CT and inflate expression
    for(gene in genes_to_sample)
    {
      # set all the expression for this celltype to 0
      counts_inflated[gene ,CT_cells] = 0
      
      gene_mean = means_perCT[grep(paste("^",gene,"$", sep=""),  rownames(means_perCT)) , which(CT == colnames(means_perCT))]
      # there are some genes that are not expressed at all in this CT -> take the mean estimated expression
      if(gene_mean < 0.001) {gene_mean = means_perCT[gene,] %>% mean}
      
      gene_dispersion = genemetadata$disp$edgeR_dispersion[which(rownames(counts) %in% gene)]
      mu = gene_mean * FC
      x1 = rnbinom(50000, mu = mu, size = 1/gene_dispersion) # shape parameter of the gamma mixing distribution
      x1 = sample(x1[x1>0] , length(cells_to_impute), replace = T)
      
      # Add the final expression to sampled zero cells
      counts_inflated[gene ,cells_to_impute] = x1
      # save the % of cells expressing the gene
      perc_cells_expressing_lst[[comb_CT]][[paste0(tmp_CT, "_" ,CT)]][[gene]] =  table(counts_inflated[gene ,CT_cells]>0)["TRUE"] / length(counts_inflated[gene ,CT_cells])
      if(is.na(perc_cells_expressing_lst[[comb_CT]][[paste0(tmp_CT, "_" ,CT)]][[gene]])) {perc_cells_expressing_lst[[comb_CT]][[paste0(tmp_CT, "_" ,CT)]][[gene]] = 0}
    
    }
  }
} 

################
# save results #
################

write.table(counts_inflated, sc_inflated_counts_path , sep = "\t")
write.table(metadata,sc_metadata_path, sep = "\t")
saveRDS(simulated_interactions_lst,simulated_interactions_path)
saveRDS(perc_cells_expressing_lst,cells_sampled_perCTCT_path)


sessionInfo()
