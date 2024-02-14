compute_diagnostic_plots = function(original_counts_aveLogCPM , master_lst, FC_param, PCE_param, dataset, CT_toPlot)
{
  plot_avelogcpm_fixed_PCE = list()
  plot_avelogcpm_fixed_FC= list()
  plot_corr_fixed_perc_cells_expressing = list()
  plot_corr_fixed_FC = list()
  
  #######################################
  ### Plot AveLogCPM having PCE fixed ###
  #######################################
  for(perc in PCE_param)
  {
    
    n = grep(paste0("^",perc,"$"), names(master_lst) %>% str_split("_") %>% lapply("[[", 4))
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
      
      current_FC = str_split( x,"_") %>% lapply(., "[[", 2) %>% unlist
      df = data.frame(original_counts = original_counts_aveLogCPM , avelogcpm = master_lst[[x]][["counts"]], is_LR =  names(master_lst[[x]][["counts"]]) %in% genes_to_plot)
      plot = ggplot(df,aes(x = original_counts , y = avelogcpm , color = is_LR)) + 
        geom_point(size = 0.5) + 
        ggtitle(paste0("PCE = " ,perc , " dataset = ",dataset, " CT = ",paste(CT_toPlot, collapse = " "))) +
        xlab("aveLogCPM original counts") +
        ylab(paste("aveLogCPM FC=",current_FC))
      plot_avelogcpm_fixed_PCE[[x]] = plot
    }
    
    ####################### % of cells expressing
    df = data.frame(FC1 = master_lst[[min_FC_index]][["PCE"]] , FC2 = master_lst[[max_FC_index]][["PCE"]] )
    plot = ggplot(df,aes(x = FC1 , y = FC2)) + 
      geom_point(size = 0.75) + 
      geom_abline(slope=1, intercept = 0) +
      ggtitle(paste("Param perc_cells_expressing is", perc, ", dataset ", dataset)) +
      xlab(paste("FC =", names(master_lst)[[min_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist)) +
      ylab(paste("FC =" , names(master_lst)[[max_FC_index]] %>% str_split("_") %>% lapply("[[",2) %>% unlist))
    plot_corr_fixed_perc_cells_expressing[[paste0("PCE_",perc)]] = plot
  }
  
  ########################################################################
  ### Plot percentage of cells expressing LR genes and having fixed FC ###
  ########################################################################
  '
  for(fc in FC_param) -> this is the same graphs as the section above
  {
    n = grep(fc, names(master_lst) %>% str_split("_") %>% lapply("[[", 2))
    # Get the filename with highest and lowers perc_cells_expressing
    min_PCE_index = n[str_split( names(master_lst[n]),"_") %>% lapply(., "[[", 4) %>% which.min]
    max_PCE_index = n[str_split( names(master_lst[n]),"_") %>% lapply(., "[[", 4) %>% which.max]
    
    ######################## Plot change in expression magnitude
    for(x in names(master_lst[n]))
    {
      current_PCE = str_split(x,"_") %>% lapply(., "[[", 4) %>% unlist
      df = data.frame(original_counts = original_counts_aveLogCPM , avelogcpm = master_lst[[x]][["counts"]], is_LR =  names(master_lst[[x]][["counts"]]) %in% genes_to_plot)
      plot = ggplot(df,aes(x = original_counts , y = avelogcpm, color = is_LR)) + 
        geom_point(size = 0.5) + 
        ggtitle(paste0("FC = " ,fc , " dataset = ",dataset, " CT = ",paste(CT_toPlot, collapse = " "))) +
        xlab("aveLogCPM original counts") +
        ylab(paste("aveLogCPM PCE=",current_PCE))
      plot_avelogcpm_fixed_FC[[x]] = plot
    }
    
    ####################### % of cells expressing
    df = data.frame(perc1 = master_lst[[min_PCE_index]][["PCE"]], perc2 = master_lst[[max_PCE_index]][["PCE"]] )
    plot = ggplot(df,aes(x = perc1 , y = perc2)) + 
      geom_point(size = 0.75) + 
      geom_abline(slope=1, intercept = 0) +
      ggtitle(paste0("FC is ",fc , ", dataset " , dataset)) +
      xlab(paste0("Param perc_cells_expressing = ", min(PCE))) +
      ylab(paste0("Param perc_cells_expressing = ", max(PCE)))
    plot_corr_fixed_FC[[paste0("FC_",fc)]] = plot
  }
  '
  return(list = list(avelogcpm_fixedPCE = plot_avelogcpm_fixed_PCE , PCE_fixedPCE = plot_corr_fixed_perc_cells_expressing))
}


plot_error = function(data, metric_plot, FC) {
  data = filter(data,metric == metric_plot)
  data$lower = data$value - data$value_sd
  data$upper = data$value + data$value_sd
  data$error_range = data$upper-data$lower 
  
  # symmetric color range
  color_range = c(min(data$error_range), max(data$error_range))
  if(FC)
  {
    p = 
      ggplot(data) +
      geom_point(aes(x=PCE,y=method, color=error_range) , size = 3) +
      # limits should be the same, using divergent palette for ease of seeing when 
      # interval contains 0
      scale_color_gradientn(colors=cetcolor::cet_pal(7, 'd1a'), limits=color_range) +
      theme_bw() +
      ggtitle(paste0("Error of ", metric_plot)) 
    
    return(p)
  } else {
    p = 
      ggplot(data) +
      geom_point(aes(x=FC,y=method, color=error_range) , size = 3) +
      # limits should be the same, using divergent palette for ease of seeing when 
      # interval contains 0
      scale_color_gradientn(colors=cetcolor::cet_pal(7, 'd1a'), limits=color_range) +
      theme_bw() +
      ggtitle(paste0("Error of ", metric_plot)) 
    
    return(p)
  }
  
}
