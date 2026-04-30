library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(ggpubr)

# load config file to determine FC and PCE parameters
config = yaml::read_yaml("config.yaml")
FC = config$semiSimulation$FC %>% unlist %>% as.double
PCE = config$semiSimulation$PCE %>% unlist %>% as.integer

# select only 3 FC and 3 PCE for the hierarchical plot
FC = c(0.3,1,10)
PCE = c(7,20,40)


# Change the file according to your dataset
master_lst = readRDS("output/results/10x_master_lst.RDS")
benchmark_simulated_genes = c(master_lst$CT1[[1]]$L , master_lst$CT2[[1]]$R) %>% unique # simulated genes are the same for entire dataset

out = compare_expressionProfile(dataset = "10x",
                                FC_values = FC,
                                PCE_values = PCE,
                                benchmark_simulated_genes = benchmark_simulated_genes)

ggarrange(plotlist = out$plots, common.legend = T)

# ----------------- #
# Functions
# ----------------- #

# compare if  distribution of specific genes in two raw count matrices are similar
# doesnt depend on the amount of cells
compare_expressionProfile = function(dataset, FC_values, PCE_values, benchmark_simulated_genes)
{
  # dataset -> dataset to retrieve counts from
  # FC_values -> FC that were used in simulation, can be retrieved from config.yaml
  # PCE_values -> PCE that were used in simulation, can be retrieved from config.yaml
  # benchmark_simulated_genes -> genes that were inflated in simulation
  
  # check for presence of columns in LR database
  if (!dataset %in% c("10x", "VASAseq", "SMARTseq2")) {
    stop(paste("Select the correct dataset variable"))
  }
  
  # Get the counts of benchmark
  file_path = file.path("output",paste0(dataset,"_semiSimulation_NB"))
  inflated_counts_files = file_path %>% list.files(., pattern = "sc_inflated_counts")
  
  
  # Get the metadata of benchmark
  metadata_files = file_path %>% list.files(., pattern = "sc_metadata")
  
  out_lst = list()
  for(PCE in PCE_values)
  {
    for(FC in FC_values)
    {
      # get the correct inflated raw counts from benchmark
      inflated_counts_file = inflated_counts_files %>% magrittr::extract(grepl(paste0("_" , FC , "_", ".*",PCE , ".tsv$"), inflated_counts_files))
      counts_benchmark = read.table(file.path(file_path,inflated_counts_file)) 
      
      # get the correct metadata
      metadata_file = metadata_files %>% magrittr::extract(grepl(paste0("_" , FC , "_", ".*",PCE , ".tsv$"), metadata_files))
      metadata = read.table(file.path(file_path,metadata_file)) 
      
      sce = SingleCellExperiment(assays = list(counts = counts_benchmark)) %>% logNormCounts # normalize
      sce = sce[, metadata %>% filter(Celltype %in% c("CT1","CT2")) %>% pull(cell_ID)]
      
      BENCH_mtx = logcounts(sce)[rownames(sce) %in% benchmark_simulated_genes,] %>% as.matrix
      
      # Convert matrix to long format
      df = melt(as.matrix(BENCH_mtx))
      colnames(df) = c("Gene", "Cell", "Expression")
      df = df %>% 
        filter(Expression > 0)
      
      # calculate quantiles
      quants = quantile(df$Expression, probs = c(0.25, 0.5, 0.75))
      
      plot = ggplot(df, aes(x = Expression)) +
        stat_ecdf(geom = "step", color = "firebrick", size = 1) +
        # Add vertical dashed lines at the calculated quantile values
        geom_vline(xintercept = quants, 
                   linetype = "dashed", color = "grey50", alpha = 0.8) +
        
        annotate("text", x = quants[2], y = 0.05, 
                 label = paste("Median:", round(quants[2], 2)), 
                 angle = 90, vjust = -0.5, size = 3, color = "grey30") +
        theme_bw() +
        labs(title = paste0("Cumulative Distribution of Expression | FC ", FC , " | PCE ", PCE),
             x = "Log-normalized expression", y = "Probability")
      
      out_lst[[paste0("PCE",PCE)]][[paste0("FC",FC)]] = plot
    }
  }
  rm(sce) ; gc()
  return(list(plots = out_lst, quantiles = quants))
}
