compute_diagnostic_plots = function(counts , master_lst, FC_param, PCE_param, dataset, cell_metadata , CT_toPlot)
{
  plot_avelogcpm_fixed_PCE = list()
  plot_avelogcpm_fixed_FC= list()
  plot_corr_fixed_PCE_cells_expressing = list()
  plot_corr_fixed_FC = list()
  
  # selects cells belonging to the celltype indicated by CT_toPlot
  original_counts = original_counts[,cell_metadata$Celltype %in% CT_toPlot]
  original_counts_aveLogCPM = aveLogCPM(original_counts)
  
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
      master_lst[[x]][["counts"]] = master_lst[[x]][["counts"]][,cell_metadata$Celltype %in% CT_toPlot]
      counts_aveLogCPM = aveLogCPM(master_lst[[x]][["counts"]])
      
      current_FC = str_split( x,"_") %>% lapply(., "[[", 2) %>% unlist
      df = data.frame(original_counts = original_counts_aveLogCPM , avelogcpm = counts_aveLogCPM, is_LR =  rownames(master_lst[[x]][["counts"]]) %in% genes_to_plot)
      plot = ggplot(df,aes(x = original_counts , y = avelogcpm , color = is_LR)) + 
        geom_point(size = 0.5) + 
        ggtitle(paste0("PCE = " ,PCE , " dataset = ",dataset, " CT = ",paste(CT_toPlot, collapse = " "))) +
        xlab("aveLogCPM original counts") +
        ylab(paste("aveLogCPM FC=",current_FC))
      plot_avelogcpm_fixed_PCE[[x]] = plot
    }
    
    ####################### % of cells expressing
    df = data.frame(FC1 = master_lst[[min_FC_index]][["PCE"]] , FC2 = master_lst[[max_FC_index]][["PCE"]] )
    plot = ggplot(df,aes(x = FC1 , y = FC2)) + 
      geom_point(size = 0.75) + 
      geom_abline(slope=1, intercept = 0) +
      ggtitle(paste("Param PCE_cells_expressing is", PCE, ", dataset ", dataset)) +
      xlab(paste("FC =", names(master_lst)[[min_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist)) +
      ylab(paste("FC =" , names(master_lst)[[max_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist))
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
