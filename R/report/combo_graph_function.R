
# Need to switch to CAS?

combo_plot_matches <- function(gd_1, gd_2,
                               thres_1, thres_2,
                               cas_key,
                               drop = TRUE,
                               gd_3=NULL,thres_3=NA,
                               grid = FALSE,
                               site_thresh = NA,
                               counts_both_sites = TRUE){
  
  guide_side_1 <- gd_1$guide_side[1]
  guide_side_2 <- gd_2$guide_side[1]
  
  gd_1 <- gd_1 %>%
    left_join(select(cas_key, CAS, chnm), by="CAS")
  
  gd_2 <- gd_2 %>%
    left_join(select(cas_key, CAS, chnm), by="CAS")
  
  if(!is.na(site_thresh)){
    chem_count <- gd_2 %>%
      select(CAS,site) %>%
      group_by(CAS) %>%
      summarise(count = n()) %>%
      filter(count >= site_thresh) %>%
      arrange(desc(count)) %>%
      pull(CAS)
    
    gd_2 <- filter(gd_2, CAS %in% chem_count)
    gd_1 <- filter(gd_1, CAS %in% chem_count)
      
  }
  
  if(all(is.null(gd_3))){
    
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
    
    if(all(levels(gd_1$Class) %in% levels(gd_2$Class))){
      
      class_order <- levels(gd_1$Class)
      
    } else {
      
      class_order <- combine_data%>%
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
        select(guide_side, CAS) %>%
        distinct() 
      
      chems_in_1_2 <- chems_in_1_2$CAS[duplicated(chems_in_1_2$CAS)]
      
      graphData_1_2 <- filter(graphData_1_2, CAS %in% chems_in_1_2)
      graphData_1_2$chnm <- factor(graphData_1_2$chnm, lev)
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
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e", "#800000", "#808000")
  
  pretty_logs_new <- toxEval:::prettyLogs(graphData_1_2$meanEAR)
 
  toxPlot_1_2 <- ggplot(data=graphData_1_2) +
    scale_y_log10(labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)  +
    geom_boxplot(aes(x=chnm, y=meanEAR, fill=Class),
                 outlier.size=1) +
    theme_bw() +
    coord_flip(clip = "off") 
  
  if(grid){
    toxPlot_1_2 <- toxPlot_1_2 +
      facet_grid(guide_up ~ guide_side, scales = "free", space = "free", labeller = label_parsed)+
      theme(strip.text.y = element_blank())
    
  } else {
    toxPlot_1_2 <- toxPlot_1_2 +
      facet_grid(. ~ guide_side, scales = "free")
  }
  
  toxPlot_1_2 <- toxPlot_1_2 +
    theme(axis.text = element_text( color = "black"),
          axis.text.y = element_text(size=9, vjust=.35),
          axis.title=element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank()) +
    guides(fill=guide_legend(ncol=6)) +
    theme(legend.position="bottom",
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line"),
          axis.ticks.y = element_blank()) +
    scale_fill_manual(values = cbValues, drop=FALSE) +
    scale_color_manual(values = c("black","red"), 
                       na.value = "black", guide = "none")
  
  if(!drop){
    toxPlot_1_2 <- toxPlot_1_2 +
      scale_x_discrete(drop=TRUE) 
  }
  
  plot_info <- ggplot_build(toxPlot_1_2)
  layout_stuff <- plot_info$layout
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    xmin_1 <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
    xmax_1 <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[2])
    xmin_2 <- 10^(layout_stuff$panel_scales_y[[2]]$range$range[1])
    xmax_2 <- 10^(layout_stuff$panel_scales_y[[2]]$range$range[2])
  } else {
    xmin_1 <- 10^(layout_stuff$panel_ranges[[1]]$x.range[1])
    xmax_1 <- 10^(layout_stuff$panel_ranges[[1]]$x.range[2])
    xmin_2 <- 10^(layout_stuff$panel_ranges[[2]]$x.range[1])
    xmax_2 <- 10^(layout_stuff$panel_ranges[[2]]$x.range[2])
  }
  
  
  if(!drop){

    if(!all(is.null(gd_3))){
      
      plot_data <- select(toxPlot_1_2$data, chnm, guide_side, guide_up) %>%
        distinct()
      
      chem_A <- unique(plot_data$chnm[plot_data$guide_up == gd_2$guide_up[1]])
      chem_B <- unique(plot_data$chnm[plot_data$guide_up == gd_3$guide_up[1]])
      
      chm_side_A <- data.frame(chnm = chem_A,
                               guide_up = gd_2$guide_up[1])
      
      chm_side_B <- data.frame(chnm = chem_B,
                               guide_up = gd_3$guide_up[1])
      
      chm_side <- bind_rows(mutate(chm_side_A, guide_side = guide_side_1),
                            mutate(chm_side_A, guide_side = guide_side_2),
                            mutate(chm_side_B, guide_side = guide_side_1),
                            mutate(chm_side_B, guide_side = guide_side_3))
      
      chm_side$chnm <- factor(chm_side$chnm, levels = levels(graphData_1_2$chnm))
      chm_side$guide_side <- factor(chm_side$guide_side, levels = levels(graphData_1_2$guide_side))

      countNonZero_1_2 <- right_join(countNonZero_1_2, chm_side, by=c("chnm","guide_side","guide_up"))
      
    } else {
      chm_side <- bind_rows(select(orderChem_1_2,-median),
                            select(orderChem_1_2,-median) ) %>%
        mutate(guide_side = factor(c(rep(guide_side_1,nrow(orderChem_1_2)),
                                     rep(guide_side_2,nrow(orderChem_1_2))),
                                   levels = c(guide_side_1, guide_side_2)))
      
      chm_side$chnm <- factor(chm_side$chnm, levels = levels(graphData_1_2$chnm))
      
      countNonZero_1_2 <- right_join(countNonZero_1_2, chm_side, by=c("chnm","Class","guide_side"))
    }


    countNonZero_1_2$nonZero[is.na(countNonZero_1_2$nonZero)] <- "*"
  }
  
  countNonZero_1_2 <- countNonZero_1_2 %>%
    mutate(ymin = ifelse(guide_side == guide_side_1, xmin_1, xmin_2),
           ymax = ifelse(guide_side == guide_side_1, xmax_1, xmax_2))

  labels_1_2 <- data.frame(x = rep(Inf,4),
                           y = c(xmin_1,thres_1,
                                 xmin_2,thres_2),
                           label = rep(c("# Sites","Threshold"),2),
                           guide_side = c(rep(guide_side_1,2),rep(guide_side_2,2)),
                           stringsAsFactors = FALSE)
  
  if(!all(is.null(gd_3))){
    labels_1_2$guide_up <- gd_2$guide_up[1] 
  }
  
  if(!all(is.na(countNonZero_1_2$hits[countNonZero_1_2$guide_side == guide_side_1]))){
    labels_1_2 <- bind_rows(labels_1_2,
                            data.frame(x = Inf, y = xmax_1, 
                                       label = "# Hits", 
                                       guide_side = guide_side_1,
                                       stringsAsFactors = FALSE))
  }
  
  if(!all(is.na(countNonZero_1_2$hits[countNonZero_1_2$guide_side == guide_side_2]))){
    labels_1_2 <- bind_rows(labels_1_2,
                            data.frame(x = Inf, y = xmax_2, 
                                       label = "# Hits", 
                                       guide_side = guide_side_2,
                                       stringsAsFactors = FALSE))
  }
  
  labels_1_2$guide_side <- factor(labels_1_2$guide_side, levels = c(guide_side_1, guide_side_2))
  
  if(!counts_both_sites){
    countNonZero_1_2 <- filter(countNonZero_1_2, guide_side == guide_side_1)
    labels_1_2 <- filter(labels_1_2, guide_side == guide_side_1)
  }
  
  toxPlot_1_2 <- toxPlot_1_2 +
    geom_text(data=countNonZero_1_2, size=2.5, 
              aes(x= chnm, label = nonZero, y=ymin)) +
    geom_text(data=countNonZero_1_2, size=2.5, 
              aes(x= chnm, label = hits, y=ymax)) +
    geom_text(data=labels_1_2, size=2.5,vjust=-0.5,
              aes(x = x,  y=y, label = label)) +
    geom_segment(data = thresh_df, aes(y = thres, yend = thres),
                 linetype="dashed", 
                 x = 1,
                 xend = Inf) 
  
  return(toxPlot_1_2)
  
}
