# Semi-simulation framework
# one can semi-simulate several combinations of CT-CT pairs
semi_simulate = function(counts , simulated_interactions_lst ,genemetadata,  metadata , combination_CT, FC, pce)
{
  perc_cells_expressing_lst = list()
  counts_inflated = counts
  not_inflated_cells = list()
  means_perCT = genemetadata$mean
  
  for(comb_CT in combination_CT)
  {
    tmp_var1 = str_split(comb_CT,"_")[[1]]
    CT_sender = tmp_var1[1]
    CT_receiver = tmp_var1[2]
    
    L_sample = simulated_interactions_lst[[comb_CT]] %>% str_split("_") %>% lapply(.,"[[",1) %>% as.character %>% setdiff("subunit") # remove subunit string from the L and R vectors
    R_sample = simulated_interactions_lst[[comb_CT]] %>% str_split("_") %>% lapply(.,"[[",2) %>% as.character %>% setdiff("subunit")
    
    # iterate over cell type combination
    for(tmp_CT in c("CTsender","CTreceiver"))  
    {
      if (tmp_CT == "CTsender" ) {genes_to_sample = L_sample ; CT = CT_sender
      } else if (tmp_CT == "CTreceiver") {genes_to_sample = R_sample ; CT = CT_receiver}
      
      # select cells belonging to CT
      CT_cells = colnames(counts)[which(metadata$Celltype == CT)]
      # only select specific percentage of cells to increase expression
      cells_to_impute = sample(CT_cells, (length(CT_cells) * pce / 100) %>% ceiling)
      #save the inflated cells for estimating mean
      not_inflated_cells[[comb_CT]][[CT]] = setdiff(CT_cells, cells_to_impute)
      # Iterate over every gene (L/R) depending on the CT and inflate expression
      for(gene_sample in genes_to_sample)
      {
        # set all the expression for this celltype to 0
        counts_inflated[gene_sample ,CT_cells] = 0
        
        gene_mean = means_perCT[grep(paste("^",gene_sample,"$", sep=""),  means_perCT$gene_names) , which(CT == colnames(means_perCT))]
        
        gene_dispersion = subset(genemetadata$disp, gene == gene_sample) %>% select(edgeR_dispersion) %>% as.numeric
        mu = gene_mean * FC
        x1 = rnbinom(1000, mu = mu, size = 1/gene_dispersion) # shape parameter of the gamma mixing distribution
        if(all(x1 == 0)) {x1 = sample(1, length(cells_to_impute), replace = T)} else {x1 = sample(x1[x1>0] , length(cells_to_impute), replace = T)}
        
        
        # Add the final expression to sampled zero cells
        counts_inflated[gene_sample ,cells_to_impute] = x1
        # save the % of cells expressing the gene
        perc_cells_expressing_lst[[comb_CT]][[paste0(tmp_CT, "_" ,CT)]][[gene_sample]] =  table(counts_inflated[gene_sample ,CT_cells]>0)["TRUE"] / length(counts_inflated[gene_sample ,CT_cells])
        if(is.na(perc_cells_expressing_lst[[comb_CT]][[paste0(tmp_CT, "_" ,CT)]][[gene_sample]])) {perc_cells_expressing_lst[[comb_CT]][[paste0(tmp_CT, "_" ,CT)]][[gene_sample]] = 0}
        
      }
    }
  } 
  
  return(list(counts_inflated = counts_inflated , perc_cells_expressing_lst = perc_cells_expressing_lst,
              not_inflated_cells = not_inflated_cells))
}

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
    
    n = grep(paste0("^",PCE,"$"), names(master_lst) %>% str_split("_") %>% lapply("[[", 4))
    # only select those files with correct FC
    n = intersect(n, which((names(master_lst) %>% str_split("_") %>% lapply("[[", 2)) %in% FC_param)) 
    # Get index of file with lowest FC
    min_FC_index = n[str_split( names(master_lst[n]),"_") %>% lapply(., "[[", 2) %>% which.min]
    max_FC_index = n[str_split( names(master_lst)[n],"_") %>% lapply(., "[[", 2) %>% which.max]
    
    ######################## Plot change in expression magnitude
    for(x in names(master_lst[n]))
    {
      # Select what genes to plot
      genes_to_plot = vector()
      for(i in CT_toPlot) {genes_to_plot = append(genes_to_plot , master_lst[[x]][[i]])}
      genes_to_plot = genes_to_plot %>% unlist %>% unique
      
      # filter the inflated counts based to contain cells belonging to the celltype indicated by CT_toPlot
      tmp_counts = master_lst[[x]][["counts"]][,metadata$Celltype %in% CT_toPlot]
      counts_aveLogCPM = aveLogCPM(tmp_counts)
      
      current_FC = str_split( x,"_") %>% lapply(., "[[", 2) %>% unlist
      df = data.frame(original_counts = filtered_original_counts_aveLogCPM , avelogcpm = counts_aveLogCPM, is_LR =  rownames(tmp_counts) %in% genes_to_plot)
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
    df = data.frame(FC1 = master_lst[[min_FC_index]][["PCE"]] , FC2 = master_lst[[max_FC_index]][["PCE"]] )
    plot = ggplot(df,aes(x = FC1 , y = FC2)) + 
      geom_point(size = 0.75) + 
      geom_abline(slope=1, intercept = 0) +
      ggtitle(paste("Param PCE_cells_expressing is", PCE, ", dataset ", dataset)) +
      xlab(paste("FC =", names(master_lst)[[min_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist)) +
      ylab(paste("FC =" , names(master_lst)[[max_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist)) +
      theme(plot.title = element_text(size=12))
    plot_corr_fixed_PCE_cells_expressing[[paste0("PCE_",PCE)]] = plot
  }
  
  return(list = list(avelogcpm_fixedPCE = plot_avelogcpm_fixed_PCE , PCE_fixedPCE = plot_corr_fixed_PCE_cells_expressing))
}

# currently not used
plot_variability = function(data, metric_plot, FC, color_range) {
  data = filter(data,metric == metric_plot)
  data$lower = data$value - data$value_sd
  data$upper = data$value + data$value_sd
  data$variability_range = data$upper-data$lower 
  
  # symmetric color range
  if(FC)
  {
    p = ggplot(data) +
      geom_point(aes(x=PCE,y=method, color=variability_range) , size = 3) +
      # limits should be the same, using divergent palette for ease of seeing when 
      # interval contains 0
      scale_color_gradientn(colors=cetcolor::cet_pal(7, 'd1a'), limits=color_range) +
      theme_bw() +
      ggtitle(paste0("variability of ", metric_plot)) 
  } else {
    p = ggplot(data) +
      geom_point(aes(x=FC,y=method, color=variability_range) , size = 3) +
      # limits should be the same, using divergent palette for ease of seeing when 
      # interval contains 0
      scale_color_gradientn(colors=cetcolor::cet_pal(7, 'd1a'), limits=color_range) +
      theme_bw() +
      ggtitle(paste0("variability of ", metric_plot)) 
  }
  return(p)
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
    ggtitle(paste0("PCE=",PCE , " | theoretical FC=",theoreticalFC , " | real FC median=",median)) +
    xlab("LR index") +
    ylab("FC after simulation (log10 scale)") +
    scale_y_log10() +
    theme(axis.text.x=element_blank(), #remove x axis labels
          plot.title = element_text(size=8)  , 
          axis.text.y = element_text(size = 8)
    )
  return(list(plot = plot, theoreticalFC = theoreticalFC, PCE_real = PCE, FC_real_median = median))
}
