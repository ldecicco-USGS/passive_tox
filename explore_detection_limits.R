loadd(tox_list)
loadd(chem_info)

tox_list$chem_data <- tox_list$chem_data %>%
  left_join(select(tox_list$chem_info,-n_sites,-Class, -chnm),
            by="CAS") %>%
  rowwise() %>%
  mutate(Value = max(`2010_MDL`, `2010_MQL`,
                     `2014_MDL`,`2014_MDL`, na.rm = TRUE)) %>%
  select(SiteID, `Sample Date`, CAS, Value) %>%
  distinct()

tox_list$chem_data <- tox_list$chem_data[is.finite(tox_list$chem_data$Value),]

cs_det_level <- get_chemical_summary(tox_list, readd(ACC), readd(filtered_ep))
cs_det_level$site <- cs_det_level$endPoint


plot_chemical_boxplots_mod <- function(chemicalSummary, 
                                       manual_remove=NULL,
                                       mean_logic = FALSE,
                                       sum_logic = TRUE,
                                       plot_ND = TRUE,
                                       font_size = NA,
                                       title = NA,
                                       pallette = NA,
                                       hit_threshold = NA){
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  chnm <- Class <- meanEAR <- x <- y <- ".dplyr"
  
  cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e", "#800000", "#808000")
  
  if(!plot_ND){
    chemicalSummary <- chemicalSummary[chemicalSummary$EAR > 0,]
  }
  
  if(length(unique(chemicalSummary$Class)) > length(cbValues)){
    n <- length(unique(chemicalSummary$Class))
    
    if(n > 20 & n<30){
      cbValues <- c(brewer.pal(n = 12, name = "Set3"),
                    brewer.pal(n = 8, name = "Set2"),
                    brewer.pal(n = max(c(3,n-20)), name = "Set1"))
    } else if (n <= 20){
      cbValues <- c(brewer.pal(n = 12, name = "Set3"),
                    brewer.pal(n =  max(c(3,n-12)), name = "Set2"))     
    } else {
      cbValues <- colorRampPalette(brewer.pal(11,"Spectral"))(n)
    }
    
  }
  
  single_site <- length(unique(chemicalSummary$site)) == 1
  
  y_label <- expression(EAR[Chem]~per~ToxCast~Assay)
  
  graphData <- toxEval:::graph_chem_data(chemical_summary = chemicalSummary, 
                                         manual_remove=manual_remove,
                                         mean_logic=mean_logic,
                                         sum_logic=sum_logic)
  
  pretty_logs_new <-  toxEval:::prettyLogs(graphData$meanEAR)
  
  countNonZero <- graphData %>%
    select(chnm, Class, meanEAR) %>%
    group_by(chnm, Class) %>%
    summarize(nonZero = as.character(sum(meanEAR>0)),
              hits = as.character(sum(meanEAR > hit_threshold)))
  
  countNonZero$hits[countNonZero$hits == "0"] <- ""
  
  label <- "# Assays"
  toxPlot_All <- ggplot(data=graphData) 
  
  if(!all(is.na(pallette))){
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=meanEAR, fill=chnm),lwd=0.1,outlier.size=1) +
      scale_fill_manual(values = pallette) +
      theme(legend.position = "none")
  } else {
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=meanEAR, fill=Class),
                   lwd=0.1,outlier.size=1)
  }
  
  toxPlot_All <- toxPlot_All +
    scale_y_log10(y_label, labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)  +
    theme_bw() +
    scale_x_discrete(drop = TRUE) +
    geom_hline(yintercept = hit_threshold, linetype="dashed", color="black") +
    theme(axis.text = element_text( color = "black"),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor =   element_blank(),
          panel.grid.major =   element_line(size =0.1),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5))  
  
  if(all(is.na(pallette))){
    toxPlot_All <- toxPlot_All +
      scale_fill_manual(values = cbValues, drop=FALSE) +
      guides(fill=guide_legend(ncol=6)) +
      theme(legend.position="bottom",
            legend.justification = "left",
            legend.background = element_rect(fill = "transparent", colour = "transparent"),
            legend.title=element_blank(),
            legend.key.height = unit(1,"line"))
  }
  
  if(!is.na(font_size)){
    toxPlot_All <- toxPlot_All +
      theme(axis.text = element_text(size = font_size),
            axis.title =   element_text(size=font_size))
  }
  
  #Saving for later!!!!
  if(packageVersion("ggplot2") >= '2.2.1.9000'){
    toxPlot_All <- toxPlot_All +
      coord_flip(clip = "off")
  } else {
    toxPlot_All <- toxPlot_All +
      coord_flip()      
  }
  
  plot_info <- ggplot_build(toxPlot_All)
  layout_stuff <- plot_info$layout
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    ymin <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
    ymax <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[2])
  } else {
    ymin <- 10^(layout_stuff$panel_ranges[[1]]$x.range[1])
    ymax <- 10^(layout_stuff$panel_ranges[[1]]$x.range[2])
  }
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=countNonZero, aes(x=chnm,label=nonZero, y=ymin), size = ifelse(is.na(font_size),2,0.30*font_size)) +
    geom_text(data=data.frame(x = Inf, y=ymin, label = label, stringsAsFactors = FALSE), 
              aes(x=x,  y=y, label = label),vjust=-0.5,
              size=ifelse(is.na(font_size),3,0.30*font_size)) 
  
  nHitsEP <- countNonZero$hits
  
  if(isTRUE(sum(as.numeric(nHitsEP), na.rm = TRUE) > 0)) {
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      geom_text(data=countNonZero, aes(x=chnm, y=ymax,label=nHitsEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
      geom_text(data=data.frame(x = Inf, y=ymax, label = "# Hits", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  if(!all(is.na(title))){
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      ggtitle(title)
    
    if(!is.na(font_size)){
      toxPlot_All_withLabels <- toxPlot_All_withLabels +
        theme(plot.title = element_text(size=font_size))
    }
  }
  
  if(!is.na(hit_threshold)) {
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      geom_text(data=data.frame(x = Inf, y=hit_threshold, label = "Threshold", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  return(toxPlot_All_withLabels)
  
}

plot_DL <- plot_chemical_boxplots_mod(cs_det_level, 
                                      font_size = 7, 
                                      title = " ")

ggsave(plot_DL, filename = "plots/detection_level_EARs.pdf", 
       width = 9, height = 11)

graphDataDL <- graph_chem_data(cs_det_level) %>% 
  arrange(desc(meanEAR))
