# Semi-simulation framework
# one can semi-simulate several combinations of CT-CT pairs
semi_simulate = function(counts , simulated_interactions_lst ,genemetadata,  metadata, FC, pce)
{
  PCE_lst = list()
  counts_inflated = counts
  not_inflated_cells = list()
  means_perCT = genemetadata$mean
  
  L_sample = simulated_interactions_lst[["ligand"]] %>% str_split("_")  %>% unlist %>% unique
  R_sample = simulated_interactions_lst[["receptor"]] %>% str_split("_")  %>% unlist %>% unique
  
  # iterate over cell type combination
  for(tmp_CT in c("CTsender","CTreceiver"))  
  {
    if (tmp_CT == "CTsender" ) {
      genes_to_sample = L_sample ; CT = "CT1"
    } else if (tmp_CT == "CTreceiver") {
      genes_to_sample = R_sample ; CT = "CT2"}
    
    # select cells belonging to CT
    CT_cells = colnames(counts)[which(metadata$Celltype == CT)]
    # only select specific percentage of cells to increase expression
    cells_to_impute = sample(CT_cells, (length(CT_cells) * pce / 100) %>% ceiling)
    #save the inflated cells for estimating mean
    not_inflated_cells[[CT]] = setdiff(CT_cells, cells_to_impute)
    # Iterate over every gene (L/R) depending on the CT and inflate expression
    for(gene_sample in genes_to_sample)
    {
      # set all the expression for this celltype to 0
      counts_inflated[gene_sample ,CT_cells] = 0
      
      gene_mean = means_perCT[grep(paste("^",gene_sample,"$", sep=""),  means_perCT$gene_names),] %>% .[CT] %>% as.numeric()
      
      gene_dispersion = genemetadata$disp %>% subset(gene == gene_sample) %>% select(edgeR_dispersion) %>% as.numeric
      mu = gene_mean * FC
      x1 = rnbinom(1000, mu = mu, size = 1/gene_dispersion) # shape parameter of the gamma mixing distribution
      # replace the expression for the CT according to the sampled values
      # in case there are not enough sampled values > 0 then sample from the > 0 values with replacement
      if(length(x1[x1>0]) > length(cells_to_impute)) {x1 = sample(x1[x1>0] , length(cells_to_impute), replace = F)
      } else if(all(x1 == 0)) {x1 = sample(1 , length(cells_to_impute), replace = T)
      } else {x1 = sample(x1[x1>0] , length(cells_to_impute), replace = T)}
      
      
      # Add the final expression to sampled zero cells
      counts_inflated[gene_sample ,cells_to_impute] = x1
      
      # save the % of cells expressing the gene
      PCE_lst[[paste0(tmp_CT, "_" ,CT)]][[gene_sample]] =  table(counts_inflated[gene_sample ,CT_cells]>0)["TRUE"] / length(counts_inflated[gene_sample ,CT_cells])
    }
  }
  
  return(list(counts_inflated = counts_inflated , PCE_lst = PCE_lst,
              not_inflated_cells = not_inflated_cells))
}

# Generate 2 plots:
# avelogcpm plot according to edgeR that shows how much signal we added to data
# Effective percentage of cells expressing the genes we simulated. This plot is just a sanity check that we are setting correctly amount of cells expressing the genes
compute_diagnostic_plots = function(counts , master_lst, FC_param, PCE_param, dataset, metadata , CT_toPlot)
{
  plot_avelogcpm_fixed_PCE = list()
  plot_avelogcpm_fixed_FC= list()
  plot_corr_fixed_PCE_cells_expressing = list()
  plot_corr_fixed_FC = list()
  
  # selects cells belonging to the celltype indicated by CT_toPlot
  filtered_raw_metadata =  filter(metadata, Celltype %in% CT_toPlot) 
  filtered_original_counts = counts[,filtered_raw_metadata$cell_ID]
  filtered_original_counts_aveLogCPM = aveLogCPM(filtered_original_counts)
  
  #######################################
  ### Plot AveLogCPM having PCE fixed ###
  #######################################
  for(PCE in PCE_param)
  {
    
    n = grep(paste0("^",PCE,"$"), names(master_lst[["PCE"]]) %>% str_split("_") %>% lapply("[[", 4))
    # only select those files with correct FC
    n = intersect(n, which((names(master_lst[["PCE"]]) %>% str_split("_") %>% lapply("[[", 2)) %in% FC_param)) 
    # Get index of file with lowest FC
    min_FC_index = n[str_split( names(master_lst[["PCE"]][n]),"_") %>% lapply(., "[[", 2) %>% which.min]
    max_FC_index = n[str_split( names(master_lst[["PCE"]])[n],"_") %>% lapply(., "[[", 2) %>% which.max]
    
    ######################## Plot change in expression magnitude
    for(x in names(master_lst[["PCE"]][n]))
    {
      # Select what genes to plot
      for(i in CT_toPlot) {genes_to_plot = master_lst[[i]][[x]] %>% unlist %>% unique}
      
      if(length(CT_toPlot) == 2) {genes_to_plot = c(master_lst$CT1[[x]][[1]], master_lst$CT2[[x]][[1]])  %>% unique}
      
      counts_aveLogCPM = master_lst[["counts_aveLogCPM"]][[x]]
      
      current_FC = str_split( x,"_") %>% lapply(., "[[", 2) %>% unlist
      df = data.frame(original_counts = filtered_original_counts_aveLogCPM , avelogcpm = counts_aveLogCPM, is_LR =  names(counts_aveLogCPM) %in% genes_to_plot)
      plot = ggplot(df,aes(x = original_counts , y = avelogcpm , color = is_LR)) + 
        geom_point(size = 0.5) + 
        ggtitle(paste0("PCE = " ,PCE , " dataset = ",dataset, " CT = ",paste(CT_toPlot, collapse = " "))) +
        xlab("aveLogCPM original counts") +
        ylab(paste("aveLogCPM FC=",current_FC)) +
        theme(plot.title = element_text(size=8) , 
              axis.text.x = element_text(size = 8) , 
              axis.text.y = element_text(size = 8))
      plot_avelogcpm_fixed_PCE[[x]] = plot
    }
    ####################### % of cells expressing
    df = data.frame(FC1 = master_lst[["PCE"]][[min_FC_index]] , FC2 = master_lst[["PCE"]][[max_FC_index]] )
    plot = ggplot(df,aes(x = FC1 , y = FC2)) + 
      geom_point(size = 0.75) + 
      geom_abline(slope=1, intercept = 0) +
      ggtitle(paste("Param PCE_cells_expressing is", PCE, ", dataset ", dataset)) +
      xlab(paste("FC =", names(master_lst[["PCE"]])[[min_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist)) +
      ylab(paste("FC =" , names(master_lst[["PCE"]])[[max_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist)) +
      theme(plot.title = element_text(size=12))
    plot_corr_fixed_PCE_cells_expressing[[paste0("PCE_",PCE)]] = plot
  }
  
  return(list = list(avelogcpm_fixedPCE = plot_avelogcpm_fixed_PCE , PCE_fixedPCE = plot_corr_fixed_PCE_cells_expressing))
}

# As theoretical FC that we apply in the semi-simulation actually doesnt represent the practical FC that the data will be transformed with generate a
# plot with real FC after semi-simulation
plot_FCafter_semisimulation = function(vec, theoreticalFC, PCE)
{
  PCE = PCE %>% unname()
  theoreticalFC = theoreticalFC %>% unname()
  df = data.frame(gene = names(vec), value = vec)
  median = median(df$value) %>% round(.,2)
  plot = ggplot(df, aes(x = gene , y = value) ) + 
    geom_point() +
    geom_hline(yintercept=theoreticalFC, linetype="dashed", color = "red", linewidth = 1)  + 
    geom_hline(yintercept=median, linetype="dashed", color = "blue", linewidth = 1)  + 
    ggtitle(paste0("PCE=",PCE , " | theoretical FC=",theoreticalFC , " | effective FC median=",median)) +
    xlab("LR index") +
    ylab("FC after simulation")  +
    coord_trans(y="log10") +
    theme(axis.text.x=element_blank(), #remove x axis labels
          plot.title = element_text(size=8)  , 
          axis.text.y = element_text(size = 8)
    )
  return(list(plot = plot, theoreticalFC = theoreticalFC, PCE_real = PCE, FC_real_median = median))
}

compare_expressionProfile = function(dataset, FC_values, PCE_values, benchmark_simulated_genes)
{
  # dataset -> dataset to retrieve counts from
  # FC_values -> FC that were used in simulation, can be retrieved from config.yaml
  # PCE_values -> PCE that were used in simulation, can be retrieved from config.yaml
  # benchmark_simulated_genes -> genes that were inflated in simulation
  
  
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
             x = "Log-normalized expression", y = "Proportion of cells in dataset")
      
      out_lst[["plots"]][[paste0("PCE",PCE)]][[paste0("FC",FC)]] = plot
      out_lst[["quantiles"]][[paste0("PCE",PCE)]][[paste0("FC",FC)]] = quants
    }
  }
  rm(sce) ; gc()
  return(out_lst)
}
