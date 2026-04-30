library(ggplot2)
library(reshape2)

# Load your count table
#counts = read.table("data/10x_immune_R1/raw_counts.tsv")
data = readRDS("/home/vkorob/Documents/data_scripts/RNA_seq/10x_scRNAseq_jessberger/Final_Annotation_sce.rds")
counts = SingleCellExperiment::counts(data)
rownames(counts) = toupper(rownames(counts))

# load config file to determine FC and PCE parameters
config = yaml::read_yaml("config.yaml")
FC = config$semiSimulation$FC %>% unlist %>% as.double
PCE = config$semiSimulation$PCE %>% unlist %>% as.integer

# Change the file according to your dataset
master_lst = readRDS("output/results/10x_master_lst.RDS")
benchmark_simulated_genes = c(master_lst$CT1[[1]]$L , master_lst$CT2[[1]]$R) %>% unique # simulated genes are the same for entire dataset

out = compare_expressionProfile(counts_USERdataset = counts,
                                dataset = "10x",
                                FC_values = FC,
                                PCE_values = PCE,
                                benchmark_simulated_genes = benchmark_simulated_genes
)

# density plots and KS statistic
ggarrange(plotlist = out$plots, common.legend = T)

# ----------------- #
# Functions
# ----------------- #

# compare if  distribution of specific genes in two raw count matrices are similar
# doesnt depend on the amount of cells
compare_expressionProfile = function(counts_USERdataset , dataset, FC_values, PCE_values, benchmark_simulated_genes ,method = "normalizedCounts", LR_database_path = "data/LR_database.tsv")
{
  # counts_USERdataset -> raw count matrix for checking
  # dataset -> dataset to retrieve counts from
  # FC_values -> FC that were used in simulation, can be retrieved from config.yaml
  # PCE_values -> PCE that were used in simulation, can be retrieved from config.yaml
  # benchmark_simulated_genes -> genes that were inflated in simulation
  # method -> normalization method. for now only 1
  # LR_database_path -> database path (.tsv) file , can be another one
  
  # check for presence of columns in LR database
  if (!dataset %in% c("10x", "VASAseq", "SMARTseq2")) {
    stop(paste("Select the correct dataset variable"))
  }
  
  # Get the counts of benchmark
  file_path = file.path("output",paste0(dataset,"_semiSimulation_NB"))
  inflated_counts_files = file_path %>% list.files(., pattern = "sc_inflated_counts")
  
  
  LRdb = read.table(LR_database_path)
  
  # check for presence of columns in LR database
  missing = setdiff(c("ligand","receptor"), colnames(LRdb))
  if (length(missing) > 0) {
    stop(paste("Error: The dataframe is missing required columns:", 
               paste(missing, collapse = ", ")))
  }
  
  # get ligand / receptor genes present in USER dataset
  LRdb_genes = LRdb %>%
    select(c(ligand,receptor)) %>% 
    unlist %>% 
    str_split("_") %>%
    unlist %>% 
    unique
  
  genes = rownames(counts_USERdataset)
  LRdb_genes = LRdb_genes[LRdb_genes %in% genes]
  
  # check if counts provided by USER has genes with lowecase (usually mouse annotations)
  check_gene_case_consistency(genes, LRdb_genes)
  
  # subset USER counts
  USER_mtx = counts_USERdataset[which(rownames(counts_USERdataset) %in% LRdb_genes), ] %>% as.matrix

  # check percentage of cells expressing ligands/receptors
  percent_expressing = mean(rowSums(USER_mtx > 0) / ncol(USER_mtx)) * 100
  
  # which value PCE from becnhmark is closest to the USER supplied counts
  PCE_use = PCE_values[which.min(abs(PCE_values - percent_expressing))]
  
  # perform KS test on distribution of non 0 values in both USER supplied and benchmark counts to find best FC from benchmark to rely on
  out_lst = list()
  for(FC in FC_values)
  {
    # get the correct inflated raw counts from benchmark
    inflated_counts_file = inflated_counts_files %>% magrittr::extract(grepl(paste0("_" , FC , "_", ".*",PCE_use , ".tsv$"), inflated_counts_files))
    counts_benchmark = read.table(file.path(file_path,inflated_counts_file))
    
    # normalize
    BENCH_mtx = counts_benchmark[benchmark_simulated_genes,] %>% as.matrix
    
    # Melt them into a 'long' format
    df1 = melt(USER_mtx)
    df2 = melt(BENCH_mtx)
    
    df1$Dataset = "USER_dataset"
    df2$Dataset = "CCC_benchmark_dataset"
    
    # combine both df and select only positive counts
    plot_df = rbind(df1, df2)
    plot_df %<>% filter(value > 0) %>%
      mutate(value = log10(value))
    
    # perform KS test to compare distributions
    statistic = ks.test(plot_df %>% filter(Dataset == "USER_dataset") %>% pull(value), plot_df %>% filter(Dataset == "CCC_benchmark_dataset") %>% pull(value))$statistic
    
    plot = ggplot(plot_df, aes(x = value, fill = Dataset)) +
      geom_density(alpha = 0.5) +
      labs(y = "density", x = "log10 counts", title = paste0("Density of > 0 counts of LR genes for FC ", FC , " and PCE ", PCE_use)) +
      theme_minimal() +
      annotate("label", x = Inf, y = Inf, 
               label = paste("KS statistic:", round(statistic, 3)),
               vjust = 1.5, hjust = 1.1, 
               fill = "white", alpha = 0.8)
    
    out_lst[[inflated_counts_file]] = plot
  }
  return(list(plots = out_lst, percent_expressing = percent_expressing))
}

# check if both gene vectors have lowercase or uppercase genes
check_gene_case_consistency = function(vec1, vec2) {
  # Determine the state of vector 1
  all_lower1 = all(vec1 == tolower(vec1))
  all_upper1 = all(vec1 == toupper(vec1))
  
  # Determine the state of vector 2
  all_lower2 = all(vec2 == tolower(vec2))
  all_upper2 = all(vec2 == toupper(vec2))
  
  # Check if both are internally consistent (not a mix of cases)
  if (!(all_lower1 | all_upper1)) stop("Vector 1 has mixed casing.")
  if (!(all_lower2 | all_upper2)) stop("Vector 2 has mixed casing.")
  
  # Check if both vectors match each other
  if (all_lower1 && all_lower2) {
    return(TRUE)
  } else if (all_upper1 && all_upper2) {
    return(TRUE)
  } else {
    # If they don't match, throw error
    stop(paste0("Casing mismatch! Vector 1 is ", 
                if(all_lower1) "lowercase" else "uppercase", 
                " but Vector 2 is ", 
                if(all_lower2) "lowercase" else "uppercase", "."))
  }
}
