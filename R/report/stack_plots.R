prep_site_list <- function(Sites){
  
  lakes_ordered <- c("Lake Superior",
                     "Lake Michigan",
                     "Lake Huron",
                     "Lake Erie",
                     "Lake Ontario")
  
  Sites$site_grouping[!(Sites$site_grouping %in% lakes_ordered)] <- "Lake Erie"
  
  Sites$site_grouping <- factor(Sites$site_grouping,
                                levels=lakes_ordered)
  
  sites_ordered <- c("StLouis","Nemadji","WhiteWI","Bad",
                     "Montreal","PresqueIsle","Pigeon","Ontonagon",
                     "Sturgeon","Tahquamenon",
                     "Manistique","Escanaba","Ford","Menominee",
                     "Peshtigo","Oconto","Fox","Manitowoc",
                     "Sheboygan #4","Sheboygan #3","Sheboygan #2","Sheboygan #1",
                     "MilwaukeeMouth","IndianaHC #1","IndianaHC #2",
                     "Burns","StJoseph","PawPaw","Kalamazoo",
                     "GrandMI #1","GrandMI #2","GrandMI #3","GrandMI #4",
                     "Muskegon","WhiteMI","PereMarquette","Manistee",
                     "Indian","Cheboygan","ThunderBay","AuSable",
                     "Rifle","Saginaw","BlackMI","Clinton",
                     "Rouge","HuronMI","Raisin",
                     "Maumee","Portage","Sandusky","HuronOH",
                     "Vermilion","BlackOH","Rocky","Cuyahoga",
                     "GrandOH","Ashtabula","Cattaraugus","Buffalo",
                     "Tonawanda","Genesee","Genesee - resampled","Oswego","BlackNY",
                     "Oswegatchie","Grass","Raquette","StRegis")
            
  Sites$`Short Name` <- factor(Sites$`Short Name`,
                               levels = sites_ordered)
  
  return(Sites)
  
}

plot_tox_stacks_manuscript <- function(chemical_summary, 
                            chem_site,
                            category = "Biological",
                            mean_logic = FALSE,
                            sum_logic = TRUE,
                            manual_remove = NULL,
                            include_legend = TRUE, 
                            font_size = NA,
                            title = NA,
                            top_num = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- maxEAR <- ".dplyr"
  SiteID <- site_grouping <- n <- index <- `Short Name` <- count <- x <- y <- label <- ".dplyr"
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- ""
  }
  
  if(category == "Chemical"){
    graphData <- graph_chem_data(chemical_summary = chemical_summary,
                                 manual_remove = manual_remove,
                                 mean_logic = mean_logic,
                                 sum_logic = sum_logic)   
    names(graphData)[names(graphData) == "maxEAR"] <- "meanEAR"
    names(graphData)[names(graphData) == "chnm"] <- "category"
  } else {
    graphData <- tox_boxplot_data(chemical_summary = chemical_summary,
                                  category = category,
                                  manual_remove = manual_remove,
                                  mean_logic = mean_logic,
                                  sum_logic = sum_logic) 
    if(category == "Chemical"){
      graphData$category <- graphData$chnm
    } 
  }

  counts <- chemical_summary %>%
    dplyr::group_by(site, CAS) %>%
    dplyr::summarize(count = sum(EAR > 0)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(count = sum(count > 0)) %>%
    dplyr::left_join(dplyr::select(chem_site, site=SiteID, `Short Name`, site_grouping), by="site") %>%
    dplyr::select(-site)
  
  siteToFind <- unique(chemical_summary$shortName)
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$category)))
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(4)
  cbValues <- sample(cbValues)
  
  siteLimits <- chem_site$`Short Name`
  single_site <- length(siteToFind) == 1
  
  if(!single_site){
    
    y_label <- toxEval:::fancyLabels(category, mean_logic, sum_logic, single_site, sep = TRUE, include_site = FALSE)
    
    graphData <- graphData %>%
      dplyr::left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                       by=c("site"="SiteID"))
    
    placement <- -0.05*diff(range(graphData$meanEAR))
    
    label_samples <- data.frame(x=-Inf,
                                y=placement,
                                label="# Detected\nChemicals", 
                                site_grouping = NA,
                                stringsAsFactors = FALSE)
    if(isTRUE(is.null(levels(chem_site$site_grouping)))){
      x <- factor(chem_site$site_grouping)
      label_samples$site_grouping <- levels(x)[1]
    } else {
      label_samples$site_grouping <- factor(levels(chem_site$site_grouping)[1],
                                            levels = levels(chem_site$site_grouping))
    }
    
    if(!is.na(top_num)){
      orig_cat <- levels(graphData$category)
      
      top_data <- graphData %>%
        dplyr::group_by(category) %>%
        dplyr::summarize(maxEAR = max(meanEAR, na.rm = TRUE)) %>%
        dplyr::arrange(dplyr::desc(maxEAR)) %>%
        dplyr::top_n(maxEAR, n=top_num) %>%
        dplyr::mutate(category = as.character(category)) %>%
        dplyr::pull(category)
      
      other_text <- paste0("Others (",length(orig_cat)-top_num,")")
      
      graphData <- graphData %>%
        dplyr::mutate(category = as.character(category),
                      category = ifelse(category %in% top_data,
                                        category, other_text),
                      category = factor(category, levels = c(top_data, other_text)))
      
    }
    
    upperPlot <- ggplot(graphData, 
                        aes(x=`Short Name`, y=meanEAR, fill = category)) +
      theme_minimal() +
      xlab("") +
      ylab(y_label[["y_label"]]) +
      facet_grid(. ~ site_grouping, scales="free", space="free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      geom_text(data = counts, 
                aes(label = count, x=`Short Name`,y = placement), 
                size=ifelse(is.na(font_size),3,0.30*font_size),inherit.aes = FALSE) +
      geom_text(data = label_samples,hjust=1,
                aes(x=x,y=y,label=label),
                size=ifelse(is.na(font_size),2,0.3*font_size),inherit.aes = FALSE) +
      labs(caption = y_label[["caption"]])  
    
  } else {
    
    y_label <- "EARs per Individual Sample"
    
    graphData <- chemical_summary %>%
      dplyr::select(-site) 
    
    placement <- -0.05*diff(range(graphData$meanEAR))
    
    dates <- dplyr::arrange(dplyr::distinct(dplyr::select(graphData, date))) 
    dates$index <- 1:(nrow(dates))
    
    graphData <- graphData %>%
      dplyr::left_join(dates, by="date")
    
    if(category == "Biological"){
      graphData$category <- graphData$Bio_category
    } else if (category == "Chemical Class"){
      graphData$category <- graphData$Class
    } else {
      graphData$category <- graphData$chnm
    }
    
    if(!is.na(top_num)){
      orig_cat <- levels(graphData$category)
      
      top_data <- graphData %>%
        dplyr::group_by(category) %>%
        dplyr::summarize(maxEAR = max(EAR, na.rm = TRUE)) %>%
        dplyr::arrange(dplyr::desc(maxEAR)) %>%
        dplyr::top_n(maxEAR, n=top_num) %>%
        dplyr::mutate(category = as.character(category)) %>%
        dplyr::pull(category)
      
      other_text <- paste0("Others (",length(orig_cat)-top_num,")")
      
      graphData <- graphData %>%
        dplyr::mutate(category = as.character(category),
                      category = ifelse(category %in% top_data,
                                        category, other_text),
                      category = factor(category, levels = c(top_data, other_text)))
      
    }
    
    upperPlot <- ggplot(graphData, aes(x=index, y=EAR, fill = category)) +
      theme_minimal() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Individual Samples") +
      ylab(y_label) 
  }
  
  upperPlot <- upperPlot +
    geom_col()  +
    theme(plot.margin = unit(c(5.5,5.5,5.5,12), "pt"))
  
  if(length(unique(graphData$category)) <= length(cbValues)){
    upperPlot <- upperPlot + 
      scale_fill_manual(name = category,values = cbValues, drop=TRUE)
    
  } 
  
  if(!include_legend){
    upperPlot <- upperPlot +
      theme(legend.position="none")
  }
  
  if(!is.na(font_size)){
    upperPlot <- upperPlot +
      theme(axis.text = element_text(size = font_size),
            strip.text = element_text(size = font_size),
            axis.title =   element_text(size=font_size))
  }
  
  if(!is.na(title)){
    upperPlot <- upperPlot +
      ggtitle(title)
    
    if(!is.na(font_size)){
      upperPlot <- upperPlot +
        theme(plot.title = element_text(size=font_size),
              plot.caption = element_text(size=font_size))
    }
  }
  
  if(utils::packageVersion("ggplot2") >= '3.0.0'){
    upperPlot <- upperPlot +
      coord_cartesian(clip = "off")
  } 
  
  return(upperPlot)
}

