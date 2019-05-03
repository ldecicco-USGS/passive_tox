combo_plot_matches_2 <- function(gd_1, gd_2,
                               thres_1, thres_2,
                               drop = TRUE,
                               gd_3=NULL,thres_3=NA,
                               grid = FALSE,
                               axis_size = 6){
  
  
  guide_side_1 <- gd_1$guide_side[1]
  guide_side_2 <- gd_2$guide_side[1]
  
  if(all(is.null(gd_3))){
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
    
    countNonZero_1_2 <- graphData_1_2 %>%
      group_by(site, chnm, Class, guide_side) %>%
      summarise(meanEAR = mean(meanEAR, na.rm=TRUE)) %>%
      group_by(chnm, Class, guide_side) %>%
      summarise(nonZero = as.character(sum(meanEAR>0)),
                hits = as.character(sum(meanEAR > ifelse(guide_side == guide_side_1, thres_1, thres_2)))) %>%
      ungroup() 

  } else {
    graphData_1_2 <- bind_rows(mutate(gd_1,
                                      Class = as.character(Class),
                                      chnm = as.character(chnm)),
                               mutate(gd_2,
                                      Class = as.character(Class),
                                      chnm = as.character(chnm)))
    
    orderChem_1_2 <- bind_rows(mutate(gd_1,
                                      Class = as.character(Class),
                                      chnm = as.character(chnm)),
                               filter(mutate(gd_2,
                                             Class = as.character(Class),
                                             chnm = as.character(chnm)),
                                      !(chnm %in% levels(gd_1$chnm)))) %>%
      group_by(CAS,Class) %>%
      summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
      ungroup() 
    
    class_order <- toxEval:::orderClass(bind_rows(gd_1, 
                                                  filter(gd_2, !(chnm %in% levels(gd_1$chnm))),
                                                  filter(gd_3, !(chnm %in% c(levels(gd_1$chnm),levels(gd_2$chnm)))))) %>%
      pull(Class)
    
    graphData_1_2 <- bind_rows(mutate(gd_1,
                                      Class = as.character(Class),
                                      chnm = as.character(chnm)),
                               mutate(gd_2,
                                      Class = as.character(Class),
                                      chnm = as.character(chnm)), 
                               mutate(gd_3,
                                      Class = as.character(Class),
                                      chnm = as.character(chnm)))
    guide_side_3 <- gd_3$guide_side[1]
    
    countNonZero_1_2 <- graphData_1_2 %>%
      group_by(site, chnm, Class, guide_side, guide_up) %>%
      summarise(meanEAR = mean(meanEAR, na.rm=TRUE)) %>%
      group_by(chnm, Class, guide_side, guide_up) %>%
      summarise(nonZero = as.character(sum(meanEAR>0)),
                hits = as.character(sum(meanEAR > ifelse(guide_side == guide_side_1, thres_1, thres_2)))) %>%
      ungroup() 
    
    class_order <- toxEval:::orderClass(bind_rows(gd_1, 
                                                  filter(gd_2, !(chnm %in% unique(gd_1$chnm))),
                                                  filter(gd_3, !(chnm %in% c(unique(gd_1$chnm),unique(gd_2$chnm))))))
  }
  
  orderChem_1_2 <- orderChem_1_2 %>%
    mutate(Class = factor(Class, levels=class_order)) %>%
    arrange(desc(Class), desc(!is.na(median)), median)
  
  graphData_1_2$Class <- factor(graphData_1_2$Class, levels = class_order)
  graphData_1_2$chnm <- factor(graphData_1_2$chnm, 
                               levels = unique(orderChem_1_2$chnm))
  
  if(drop){
    if(grid){
      chems_in_A <- filter(graphData_1_2, guide_up == gd_2$guide_up[1]) %>%
        select(guide_side, chnm) %>%
        distinct() 
      
      chems_in_A <- chems_in_A$chnm[duplicated(chems_in_A$chnm)]
      if(length(chems_in_B) > 0){
        chems_in_A <- data.frame(chnm = chems_in_A, guide_up = gd_2$guide_up[1], stringsAsFactors = FALSE)
      } else {
        chems_in_A <- data.frame(chnm = NA, guide_up = gd_2$guide_up[1], stringsAsFactors = FALSE)
      }
      
      chems_in_B <- filter(graphData_1_2, guide_up == gd_3$guide_up[1]) %>%
        select(guide_side, chnm) %>%
        distinct()
      
      chems_in_B <- chems_in_B$chnm[duplicated(chems_in_B$chnm)]
      if(length(chems_in_B) > 0){
        chems_in_B <- data.frame(chnm=chems_in_B, guide_up = gd_3$guide_up[1], stringsAsFactors = FALSE)
      } else {
        chems_in_B <- data.frame(chnm = NA, guide_up = gd_2$guide_up[1], stringsAsFactors = FALSE)
      }
      
      chems_in <- bind_rows(chems_in_A, chems_in_B)
      
      graphData_1_2 <- graphData_1_2 %>%
        right_join(chems_in, by=c("chnm","guide_up"))
    } else {
      chems_in_1_2 <- graphData_1_2 %>%
        select(guide_side, chnm) %>%
        distinct() 
      
      chems_in_1_2 <- chems_in_1_2$chnm[duplicated(chems_in_1_2$chnm)]
      
      graphData_1_2 <- filter(graphData_1_2, chnm %in% chems_in_1_2$chnm)
    }
    
  }

  thresh_df <- data.frame(guide_side = c(guide_side_1,guide_side_2),
                          thres = c(thres_1, thres_2),
                          stringsAsFactors = FALSE)
  
  thresh_df <- thresh_df[!is.na(thresh_df$thres),]
  
  graphData_1_2$guide_side <- factor(graphData_1_2$guide_side, levels = c(guide_side_1,guide_side_2))
  thresh_df$guide_side <- factor(thresh_df$guide_side, levels = c(guide_side_1,guide_side_2))
  countNonZero_1_2$guide_side <- factor(countNonZero_1_2$guide_side, levels = c(guide_side_1,guide_side_2))
  
  cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#003366","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e", "#800000", "#808000")
  
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
    guides(fill=guide_legend(ncol=5)) +
    theme(legend.position="bottom",
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=4),
          legend.key.height = unit(1,"line"),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_fill_manual(values = cbValues, drop=FALSE) +
    scale_color_manual(values = c("black","red"), 
                       na.value = "black", guide = "none")
  if(!drop){
    toxPlot_1_2 <- toxPlot_1_2 +
      scale_x_discrete(drop=TRUE) 
  }
  
  if(grid){
    toxPlot_1_2 <- toxPlot_1_2 +
      facet_grid(guide_up ~ guide_side, 
                 scales = "free",
                 space = "free_y", 
                 labeller = label_parsed)+
      theme(strip.text.y = element_blank())
    
  } else {
    toxPlot_1_2 <- toxPlot_1_2 +
      facet_grid(. ~ guide_side, scales = "free")
  }
  
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

add_label <- function(plot_full, labels_to_use = c("A","B")){
  
  plot_data <- plot_full$data

  textData <- plot_data %>%
    group_by(guide_up, guide_side) %>%
    summarize(chnm = tail(levels(chnm),1),
              y = max(meanEAR, na.rm = TRUE)) %>%
    mutate(textExplain = labels_to_use)
  
  plot_full_w_label <- plot_full +
    geom_text(data = textData, 
              aes(label = textExplain, x = chnm, y=y),size = 2)
  
  return(plot_full_w_label)
}

strip_graph <- function(plot_full){
  
  no_axis <- plot_full +
    theme(axis.title.y=element_blank(),
          strip.text.x = element_text(size = 5),
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
    right_join(distinct(select(plot_data, CAS, Class, chnm, guide_up)), by="CAS") %>%
    mutate(guide_side = "Sites")
  
  site_counts$nSites[is.na(site_counts$nSites)] <- 0
  
  return(site_counts)
}

site_count_plot <- function(site_counts, axis_size = 6){
  
  pretty_logs_new <- toxEval:::prettyLogs(c(10,100,1000))
  
  site_graph <- ggplot(data = site_counts) +
    geom_text(aes(x=chnm, label = nSites, y=100), size = 2) +
    facet_grid(guide_up ~ guide_side, scales = "free", space = "free_y")+
    theme_bw() +
    coord_flip() +
    scale_y_log10(labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)  +
    theme(axis.text.x = element_text(size=axis_size, color = "transparent"),
          axis.text.y = element_text(size=axis_size, vjust=.35, color = "black"),
          axis.title=element_blank(),
          strip.text.x = element_text(size = 5),
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

fancy_combo <- function(graphData_tox, graphData_conc, tox_list, axis_size = 6){
  
  toxPlot_wq <- combo_plot_matches_2(graphData_tox, 
                                     graphData_conc, 
                                     thres_1 = NA, 
                                     thres_2 = NA, 
                                     drop = FALSE,
                                     axis_size = axis_size)
  
  toxPlot_wq_w_lab <- add_label(toxPlot_wq)
  no_axis <- strip_graph(toxPlot_wq_w_lab)
  site_counts_df <- site_counts(tox_list$chem_data, no_axis$data)
  site_graph <- site_count_plot(site_counts_df)
    
  return(list(no_axis=no_axis, site_graph=site_graph))
}

