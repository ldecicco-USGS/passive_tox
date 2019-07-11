
plot_tox_endpoints_manuscript <- function(chemicalSummary, AOP,
                              category = "Biological",
                              filterBy = "All",
                              manual_remove = NULL,
                              hit_threshold = NA,
                              mean_logic = FALSE, 
                              sum_logic = TRUE,
                              font_size = NA,
                              title = NA,
                              pallette = NA,
                              threshold = 0.001,
                              siteThreshold = 10){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  site <- endPoint <- EAR <- sumEAR <- meanEAR <- x <- y <- ".dplyr"
  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else if(category == "Chemical Class") {
    chemicalSummary$category <- chemicalSummary$Class
  } else {
    chemicalSummary$category <- chemicalSummary$chnm
  }
  
  single_site <- length(unique(chemicalSummary$site)) == 1
  
  if(filterBy != "All"){
    if(!(filterBy %in% unique(chemicalSummary$category))){
      stop("filterBy argument doesn't match data")
    }
    
    chemicalSummary <- chemicalSummary %>%
      filter_(paste0("category == '", filterBy,"'"))
  }
  
  endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
    group_by(endPoint,site,date) %>%
    summarize(EARsum = sum(EAR)) %>%
    group_by(site,endPoint) %>%
    summarize(EARmax = max(EARsum)) %>%
    filter(EARmax >= threshold) %>%
    group_by(endPoint) %>%
    summarize(numSites = n_distinct(site)) %>%
    arrange(desc(numSites)) %>%
    filter(numSites >= siteThreshold)
  
  priority_endpoints <- endpoints_sites_hits$endPoint
  
  chemicalSummaryPriority <- filter(chemicalSummary, endPoint %in% priority_endpoints)
  
  eps_with_ids <- unique(AOP$endPoint)
  
  chemicalSummaryPriority$has_AOP <- "AOP Undefined"
  chemicalSummaryPriority$has_AOP[chemicalSummaryPriority$endPoint %in% eps_with_ids] <- "AOP Associated"
  
  y_label <- toxEval:::fancyLabels(category, mean_logic, sum_logic, single_site)
  
  graphData <- chemicalSummaryPriority %>%
    group_by(site,date,endPoint,has_AOP) %>%
    summarise(sumEAR=sum(EAR)) %>%
    group_by(site, endPoint,has_AOP) %>%
    summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) 
    
  pretty_logs_new <- toxEval:::prettyLogs(graphData$meanEAR)
    
  chem_ns <- chemicalSummaryPriority %>%
    group_by(endPoint) %>%
    summarize(nChems = length(unique(chnm)))

  countNonZero <- graphData %>%
    group_by(endPoint) %>%
    summarise(nonZero = as.character(sum(meanEAR>0)),
              hits = as.character(sum(meanEAR > hit_threshold)))

  orderColsBy <- graphData %>%
    group_by(endPoint) %>%
    summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
    arrange(median)

  orderedLevelsEP <- orderColsBy$endPoint

  if(any(is.na(orderColsBy$median))){
    orderedLevelsEP <- c(orderColsBy$endPoint[is.na(orderColsBy$median)],
                        orderColsBy$endPoint[!is.na(orderColsBy$median)])
  }

  
  graphData$column <- ""
  
  counts_df <- countNonZero %>%
    select(-hits) %>%
    mutate(site = NA,
           has_AOP = NA,
           meanEAR = NA,
           column = "Sites")
  
 nChem_df <- chem_ns %>%
    mutate(site = NA,
           has_AOP = NA,
           meanEAR = NA,
           column = "Chemicals")
 
  AOP_IDs <- AOP %>%
    filter(endPoint %in% nChem_df$endPoint) %>%
    group_by(endPoint) %>%
    summarize(AOP_ID = length(unique(ID))) %>%
    mutate(column = "AOPs",
           site = NA,
           has_AOP = NA,
           meanEAR = NA)
  
  graphData <- bind_rows(graphData, counts_df, nChem_df, AOP_IDs)
  
  graphData$endPoint <- factor(graphData$endPoint, levels = orderedLevelsEP)
  graphData$column <- factor(graphData$column, 
                             levels = c("Sites",
                                        "",
                                        "Chemicals",
                                        "AOPs"))
  
  stackedPlot <- ggplot()+
    scale_y_log10(expression(EAR[mixture]),labels=toxEval:::fancyNumbers,breaks=pretty_logs_new) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.1),
          axis.title = element_blank()) +
    geom_boxplot(aes(x=endPoint, y=meanEAR, fill = has_AOP),
                 data = filter(graphData, column == ""), 
                 outlier.size=0.25, lwd=0.01, fatten=1) +
    coord_flip()
  
  if(!all(is.na(pallette))){
    stackedPlot <- stackedPlot +
      scale_fill_manual(values = pallette) 
  }
  
  if(is.na(title)){
    title <- "  "
  } 
  
  stackedPlot <- stackedPlot +
      ggtitle(title)    
  
  if(!is.na(font_size)){
    stackedPlot <- stackedPlot +
      theme(plot.title = element_text(size=font_size))
  }
  
  count_plot <- ggplot() +
    geom_text(aes(x=endPoint, y=1, label = nonZero),size=5*font_size/14,
                 data = filter(graphData, column == "Sites")) +
    theme_minimal() +
    theme(axis.text.y = element_text(vjust = .25,hjust=1),
          legend.position = "none",
          axis.text.x = element_text(color = "white"),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = font_size)) +
    coord_flip() +
    ggtitle("Sites")
  
  chem_plot <- ggplot() +
    geom_text(aes(x=endPoint, y=1, label = nChems),size=5*font_size/14,
            data = filter(graphData, column == "Chemicals")) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(color = "white"),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = font_size)) +
    coord_flip() +
    ggtitle("Chemicals")
  
  aop_plot <- ggplot() +
    geom_text(aes(x=endPoint, y=1, label = AOP_ID),size=5*font_size/14,
              data = filter(graphData, column == "AOPs")) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(color = "white"),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = font_size)) +
    coord_flip() +
    scale_x_discrete(drop=FALSE) +
    ggtitle("AOPs")
  
  if(!is.na(font_size)){
    stackedPlot <- stackedPlot +
      theme(axis.text = element_text(size = font_size))
    
    count_plot <- count_plot +
      theme(axis.text = element_text(size = font_size))
  }
  
  return(list(count_plot=count_plot,
              stackedPlot=stackedPlot,
              chem_plot=chem_plot,
              aop_plot=aop_plot))
}