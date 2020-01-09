combine_gd <- function(gd_1, gd_2){
  
  guide_side_1 <- gd_1$guide_side[1]
  guide_side_2 <- gd_2$guide_side[1]
  
  if(length(levels(gd_1$chnm)) > 0){
    orderChem_1_2 <- bind_rows(gd_1,
                               filter(gd_2, !(chnm %in% levels(gd_1$chnm)))) %>%
      group_by(chnm,Class) %>%
      summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
      ungroup()
  } else {
    combine_data <- gd_2 %>%
      filter(!(CAS %in% unique(gd_1$CAS))) %>%
      mutate(Class = as.character(Class),
             chnm = as.character(chnm)) %>%
      bind_rows(mutate(gd_1, 
                       Class = as.character(Class),
                       chnm = as.character(chnm)))
    
    orderChem_1_2 <- combine_data %>%
      group_by(chnm,Class) %>%
      summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
      ungroup()
  }
  
  if(all(levels(gd_2$Class) %in% levels(gd_1$Class)) &
     all(levels(gd_1$Class) %in% levels(gd_2$Class))){
    class_order <- levels(gd_1$Class)
  } else {
    class_order <- combine_data %>%
      toxEval:::orderClass() %>%
      pull(Class)
  }
  
  graphData_1_2 <- bind_rows(mutate(gd_1,
                                    Class = as.character(Class),
                                    chnm = as.character(chnm)),
                             mutate(gd_2,
                                    Class = as.character(Class),
                                    chnm = as.character(chnm)))
  
  orderChem_1_2 <- orderChem_1_2 %>%
    mutate(Class = factor(Class, levels=class_order)) %>%
    arrange(desc(Class), desc(!is.na(median)), median)
  
  graphData_1_2$Class <- factor(graphData_1_2$Class, levels = class_order)
  graphData_1_2$chnm <- factor(graphData_1_2$chnm, 
                               levels = unique(orderChem_1_2$chnm))
  
  graphData_1_2$guide_side <- factor(graphData_1_2$guide_side, levels = c(guide_side_1,guide_side_2))
  
  return(graphData_1_2)
}

class_colors <- function(chemicalSummary){
  
  # classes <- levels(chemicalSummary$Class)
  
  classes <- c("Insecticide","Flavor/Fragrance",
              "Antimicrobial disinfectant", "Herbicide",
              "Fire retardant", "Detergent metabolites",
              "Pharmaceuticals", "Plasticizer",
              "WW", "PAHs",
              "Other", "OC Pesticides",
              "Food Additive/Plasticizer", "Dye/Pigment",
              "Solvent", "PBDEs", "Sterol", "Fuel", "PCBs")

  
  cbValues <- c("brown1", "gold","darkred","darkblue","yellow","grey10",
                "darkolivegreen","darksalmon", "darkolivegreen1","cyan3","deeppink2","grey50",
                "aquamarine","azure2","darkorange","darkorchid","cornflowerblue","cornsilk","white")

  names(cbValues) <- classes
  

  return(cbValues)
  
}

combo_plot_matches_2 <- function(graphData_1_2,
                               axis_size = 6,
                               cbValues){
  
  pretty_logs_new <- toxEval:::prettyLogs(graphData_1_2$meanEAR)
  
  toxPlot_1_2 <- ggplot(data=graphData_1_2) +
    scale_y_log10(labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)  +
    geom_boxplot(aes(x=chnm, y=meanEAR, fill=Class),
                 outlier.size=0.25, lwd=0.01, fatten=1) +
    theme_bw() +
    coord_flip() +
    theme(axis.text = element_text(size=axis_size, color = "black"),
          axis.text.y = element_text(vjust=.35),
          axis.title=element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 0.05)) +
    guides(fill=guide_legend(ncol=3)) +
    theme(legend.position="bottom",
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=5),
          legend.key.height = unit(1,"line"),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.1)) +
    scale_fill_manual(values = cbValues, drop=FALSE) +
    scale_color_manual(values = c("black","red"), 
                       na.value = "black", guide = "none") +
    facet_grid(. ~ guide_side, scales = "free") 

  return(toxPlot_1_2)
  
}

graph_chem_data_CAS <- function(chemical_summary, 
                                manual_remove=NULL,
                                mean_logic = FALSE,
                                sum_logic = TRUE){
  
  site <- chnm <- Class <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  chemical_summary <- chemical_summary %>%
    select(-chnm) %>%
    distinct()
  
  if(!sum_logic){
    graphData <- chemical_summary %>%
      dplyr::group_by(site, CAS, Class) %>%
      dplyr::summarise(meanEAR=ifelse(mean_logic,mean(EAR),max(EAR))) %>%
      dplyr::ungroup()     
  } else {
    #With new dplyr...will need to filter out na's in meanEAR
    graphData <- chemical_summary %>%
      dplyr::group_by(site,date,CAS,Class) %>%
      dplyr::summarise(sumEAR=sum(EAR,na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(site, CAS,Class) %>%
      dplyr::summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
      dplyr::ungroup() 
  }
  
  if(!is.null(manual_remove)){
    graphData <- dplyr::filter(graphData, !(chnm %in% manual_remove))
  }
  
  return(graphData)
}

add_astrict <- function(plot_full, graphData_1_2){
  
  countNonZero_1_2 <- graphData_1_2 %>%
    group_by(site, chnm, Class, guide_side) %>%
    summarise(meanEAR = mean(meanEAR, na.rm=TRUE)) %>%
    group_by(chnm, Class, guide_side) %>%
    summarise(nonZero = as.character(sum(meanEAR>0))) %>%
    ungroup() 

  min_ys <- graphData_1_2 %>%
    filter(meanEAR > 0) %>%
    group_by(guide_side) %>%
    summarize(y = min(meanEAR, na.rm = TRUE))
  
  chm_side_A <- graphData_1_2 %>%
    select(chnm, Class) %>%
    distinct() %>%
    mutate(guide_side = levels(graphData_1_2$guide_side)[1],
           guide_side = factor(guide_side, levels = levels(graphData_1_2$guide_side)))

  chm_side_B <- graphData_1_2 %>%
    select(chnm, Class) %>%
    distinct() %>%
    mutate(guide_side = levels(graphData_1_2$guide_side)[2],
           guide_side = factor(guide_side, levels = levels(graphData_1_2$guide_side)))
  
  chm_side <- bind_rows(chm_side_A, chm_side_B) %>%
    left_join(min_ys, by="guide_side")
  
  countNonZero_1_2 <- right_join(countNonZero_1_2, 
                                 chm_side, 
                                 by=c("chnm","Class","guide_side"))
  
  countNonZero_1_2$nonZero[is.na(countNonZero_1_2$nonZero)] <- "*"
  
  countNonZero_1_2 <- countNonZero_1_2 %>%
    filter(nonZero == "*")
  
  plot_full_w_ast <- plot_full +
    geom_text(data=countNonZero_1_2, size=2.5, color="blue",
              aes(x= chnm, label = nonZero, y = y)) 
  
  return(plot_full_w_ast)
  
}

label_info <- function(gd1_2, labels_to_use = c("A","B")){
  
  textData <- gd1_2 %>%
    group_by(guide_side) %>%
    summarize(chnm = factor(tail(levels(chnm),1), levels = levels(chnm)),
              y = max(meanEAR, na.rm = TRUE)) %>%
    mutate(textExplain = labels_to_use)
  
  return(textData)
}

add_label <- function(plot_full, label_info, label_size = 2){
  
  plot_full_w_label <- plot_full +
    geom_text(data = label_info, 
              aes(label = textExplain, x = chnm, y=y),size = label_size)
  
  return(plot_full_w_label)
}

strip_graph <- function(plot_full, font_size = 5){
  
  no_axis <- plot_full +
    theme(axis.title.y=element_blank(),
          strip.text.x = element_text(size = font_size),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0.25, 0, 0), "cm"))
  
  return(no_axis)
}

site_counts <- function(chem_data, plot_data){
  
  site_counts <- chem_data %>%
    filter(Value > 0) %>%
    group_by(CAS) %>%
    summarise(nSites = length(unique(SiteID))) %>%
    ungroup() %>%
    right_join(distinct(select(plot_data, CAS, Class, chnm)), by="CAS") %>%
    mutate(guide_side = "Sites")
  
  site_counts$nSites[is.na(site_counts$nSites)] <- 0
  
  return(site_counts)
}

site_count_plot <- function(site_counts, axis_size = 6){
  
  pretty_logs_new <- toxEval:::prettyLogs(c(10,100,1000))
  
  site_graph <- ggplot(data = site_counts) +
    geom_text(aes(x=chnm, label = nSites, y=100), size = axis_size*5/14) +
    facet_grid(. ~ guide_side, scales = "free", space = "free_y")+
    theme_bw() +
    coord_flip() +
    scale_y_log10(labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)  +
    theme(axis.text.x = element_text(size=axis_size, color = "transparent"),
          axis.text.y = element_text(size=axis_size, vjust=.35, color = "black"),
          axis.title=element_blank(),
          strip.text.x = element_text(size = axis_size),
          plot.margin = unit(c(0, 0, 0, 0.25), "cm"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(fill = "transparent", color = "transparent"),
          axis.ticks.x = element_line(color = "transparent")) 
  
  return(site_graph)
}

fancy_combo <- function(gd1, gd2, tox_list, axis_size = 6){
  
  gd1_2 <- combine_gd(gd_1 = gd1, 
                      gd_2 = gd2)
  
  color_map <- class_colors(tox_list)
  
  toxPlot_wq <- combo_plot_matches_2(gd1_2,
                                     axis_size = axis_size,
                                     color_map)
  ast_plot <- add_astrict(toxPlot_wq, gd1_2)
  text_df <- label_info(gd1_2)
  toxPlot_wq_w_lab <- add_label(ast_plot, text_df)
  no_axis <- strip_graph(toxPlot_wq_w_lab)
  site_counts_df <- site_counts(tox_list$chem_data, no_axis$data)
  site_graph <- site_count_plot(site_counts_df, axis_size = axis_size)
    
  return(list(no_axis=no_axis, 
              site_graph=site_graph))
}

