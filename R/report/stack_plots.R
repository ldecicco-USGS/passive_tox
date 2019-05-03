prep_site_list <- function(Sites){
  
  lakes_ordered <- c("Lake Superior",
                     "Lake Michigan",
                     "Lake Huron",
                     "Lake Erie",
                     "Lake Ontario")
  
  Sites$site_grouping[!(Sites$site_grouping %in% lakes_ordered)] <- "Lake Erie"
  
  Sites$site_grouping <- factor(Sites$site_grouping,
                                levels=lakes_ordered)
  Sites$`Short Name`[Sites$`Short Name` == "Genesee - resampled"] <- "GeneseeDock"
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
                     "Tonawanda","Genesee","GeneseeDock","Oswego","BlackNY",
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
    dplyr::select(-site) %>%
    dplyr::mutate(count= ifelse(nchar(count) < 2, paste0(" ",count), count))
  
  siteToFind <- unique(chemical_summary$shortName)
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$category)))
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(4)
  cbValues <- sample(cbValues)
  
  siteLimits <- chem_site$`Short Name`
  single_site <- length(siteToFind) == 1
  
  y_label <- toxEval:::fancyLabels(category, mean_logic, sum_logic, single_site, sep = TRUE, include_site = FALSE)
  
  graphData <- graphData %>%
    dplyr::left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                     by=c("site"="SiteID"))
  graphData$`Short Name` <- factor(graphData$`Short Name`, levels = rev(levels(graphData$`Short Name`)))
  counts$`Short Name` <- factor(counts$`Short Name`, levels = rev(levels(counts$`Short Name`)))
  
  placement <- -0.05*diff(range(graphData$meanEAR))
  
  label_samples <- data.frame(x=Inf,
                              y=-Inf,
                              label="# Chemicals", 
                              site_grouping = NA,
                              stringsAsFactors = FALSE)
  
  label_samples$site_grouping <- factor(levels(chem_site$site_grouping)[1],
                                        levels = levels(chem_site$site_grouping))
  
  
  upperPlot <- ggplot(graphData) +
    geom_col(aes(x=`Short Name`, y=meanEAR, fill = category))  +
    theme_minimal() +
    xlab("") +
    # ylab(y_label[["y_label"]]) +
    facet_grid(site_grouping ~ ., scales="free", space="free") +
    geom_text(data = counts, hjust = 0,
              aes(label = count, x=`Short Name`),y = -Inf, 
              size=ifelse(is.na(font_size),2.5,0.28*font_size)) +
    geom_text(data = distinct(select(counts, site_grouping)),
              aes(label = site_grouping),
              size=ifelse(is.na(font_size),5,0.3*font_size),
              x = Inf, y = Inf, hjust = 1, vjust = 1.5) +
    geom_text(data = label_samples,
              aes(x=x,y=y,label=label),hjust=0,vjust=0,
              size=ifelse(is.na(font_size),1.75,0.3*font_size)) +
    # labs(caption = y_label[["caption"]]) +
    coord_flip(clip = "off") + 
    scale_fill_manual(name = category,
                      values = cbValues, drop=TRUE) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x = element_blank(),
          legend.position="bottom",
          legend.text = element_text(size = ifelse(is.na(font_size),2.5,0.8*font_size)),
          legend.title=element_blank(),
          legend.key.size = unit(0.35, "cm"),
          legend.box.margin=margin(-15,0,0,0)) +
    guides(fill=guide_legend(ncol=4)) 

  upperPlot
  
  if(!is.na(font_size)){
    upperPlot <- upperPlot +
      theme(axis.text = element_text(size = font_size),
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

  
  return(upperPlot)
}

