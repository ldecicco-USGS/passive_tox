chem_counts <- function(chemical_summary, chem_site){
  
  counts <- chemical_summary %>%
    dplyr::group_by(site, CAS) %>%
    dplyr::summarize(count = sum(EAR > 0)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(count = sum(count > 0)) %>%
    dplyr::full_join(dplyr::select(chem_site, site=SiteID, `Short Name`, site_grouping), by="site") %>%
    dplyr::select(-site) %>%
    dplyr::mutate(count= ifelse(nchar(count) < 2, paste0(" ",count), count),
                  count_title = "Chemicals")
  
  counts$`Short Name` <- factor(counts$`Short Name`, levels = rev(levels(counts$`Short Name`)))
  
  return(counts)
  
}

chem_count_plot <- function(chem_df, axis_size){
  
  pretty_logs_new <- toxEval:::prettyLogs(c(10,100,1000))
  
  upperPlot <- ggplot(chem_df) +
    geom_text(aes(x=`Short Name`, label = count, y=100), size = 5*axis_size/14)  +
    theme_minimal() +
    ylab("") +
    facet_grid(site_grouping ~ count_title, scales="free", space="free") +
    coord_flip(clip = "off") + 
    theme(strip.background = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_text(size = 5),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.title.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0.25), "cm"),
          axis.text.x = element_text(size = axis_size, color = "transparent"),
          axis.text.y = element_text(size = axis_size),
          axis.title.x = element_text(size = axis_size))

  return(upperPlot)
  
}

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
                            chem_site,cbValues,
                            category = "Biological",
                            mean_logic = FALSE,
                            sum_logic = TRUE,
                            manual_remove = NULL,
                            include_legend = TRUE, 
                            font_size = NA,
                            title = NA){
  
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

  graphData <- graphData %>%
    dplyr::full_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                     by=c("site"="SiteID"))
  
  graphData$`Short Name` <- factor(graphData$`Short Name`, levels = rev(levels(graphData$`Short Name`)))
  
  if(is.na(title)){
    graphData$count_title <- ""
  } else {
    graphData$count_title <- title
  }

  upperPlot <- ggplot(graphData) +
    geom_col(aes(x=`Short Name`, y=meanEAR, fill = category))  +
    theme_minimal() +
    # ylab("Sum of Maximum EAR Values") +
    facet_grid(site_grouping ~ count_title, scales="free", space="free") +
    # geom_text(data = distinct(select(graphData, site_grouping)),
    #           aes(label = site_grouping),
    #           size=5*font_size/14,
    #           x = Inf, y = Inf, hjust = 1, vjust = 1.5) +
    coord_flip(clip = "off") + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    scale_fill_manual(name = category,
                      values = cbValues, drop=TRUE) +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = font_size),
          axis.title.y = element_blank(),
          legend.position="bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.1),
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=5),
          legend.key.height = unit(0.5,"line"),
          legend.key.width = unit(0.5, "line")) +
    # guides(fill=guide_legend(ncol=5)) +
    theme(axis.text = element_text(size = font_size),
          axis.title =   element_text(size= font_size),
          axis.ticks = element_line(colour = "black", size = 0.05))

    return(upperPlot)
}

whole_stack <- function(chemicalSummary, 
                        site_info, tox_list, 
                        color_map, 
                        font_size, 
                        category, 
                        title){
  
  plot_back <- plot_tox_stacks_manuscript(chemicalSummary,
                                          cbValues = color_map, 
                                          chem_site = site_info,
                                          font_size = font_size,
                                          category = category,
                                          title = title)
  no_axis_plot_back <- strip_graph(plot_back, font_size)
  chem_df <- chem_counts(chemicalSummary, site_info)
  chem_count_graph <- chem_count_plot(chem_df,
                                      axis_size = font_size)
  
  return(list(chem_count=chem_count_graph,no_axis=no_axis_plot_back))
}


plot_tox_stacks_manuscript2 <- function(chemical_summary, 
                                       chem_site,cbValues,
                                       category = "Biological",
                                       mean_logic = FALSE,
                                       sum_logic = TRUE,
                                       manual_remove = NULL,
                                       include_legend = TRUE, 
                                       font_size = NA,
                                       title = NA){
  
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
  
  graphData <- graphData %>%
    dplyr::full_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                     by=c("site"="SiteID"))
  
  graphData$`Short Name` <- factor(graphData$`Short Name`, levels = rev(levels(graphData$`Short Name`)))
  
  if(is.na(title)){
    graphData$count_title <- ""
  } else {
    graphData$count_title <- title
  }
  
  counts_df <- chem_counts(chemical_summary, chem_site)
  
  counts_df <- counts_df %>% 
    right_join(select(chem_site, `Short Name`, map_nm), by="Short Name")
  
  counts_df <- counts_df[!duplicated(counts_df$`Short Name`),]
  
  labels_df <- data.frame(y = c(-0.05,-0.01),
                          x = c(Inf,Inf),
                          label = c("Map Name","Chemicals"),
                          site_grouping = c("Lake Superior","Lake Superior"))
  
  classes <- c("Insecticide","Flavor/Fragrance",
               "Antimicrobial disinfectant", "Herbicide",                 
               "Fire retardant", "Detergent metabolites",  
               "Pharmaceuticals", "Plasticizer",               
               "WW", "PAHs",                      
               "Other", "OC Pesticides",             
               "Food Additive/Plasticizer", "Dye/Pigment",             
               "Solvent", "PBDEs", "Sterol", "Fuel")
  
  graphData$category <- factor(graphData$category,
                               levels = classes)
  
  upperPlot <- ggplot() +
    geom_col(data = graphData, 
             aes(x=`Short Name`, y=meanEAR, fill = category))  +
    theme_minimal() +
    ylab("Sum of Maximum EAR Values") +
    geom_text(data = counts_df, 
              aes(x=`Short Name`, label = count, y=-0.01), 
              hjust = 0.5, vjust = 0.35, size = font_size/2.5) +
    geom_text(data = counts_df, 
              aes(x=`Short Name`, label = map_nm, y=-0.05), 
              hjust = 0.5, vjust = 0.35, size = font_size/2.5) +
    geom_text(data = labels_df,
              aes(y = y, x = x, label = label), 
              vjust = -0.5, size = font_size/3) +
    facet_grid(site_grouping ~ ., scales="free", space="free") +
    coord_flip(clip = "off") + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    scale_fill_manual(name = category,
                      values = cbValues, drop=TRUE) +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = font_size),
          axis.title.y = element_blank(),
          legend.position="bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.1),
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=5),
          legend.key.height = unit(0.5,"line"),
          legend.key.width = unit(0.5, "line"),
          axis.text = element_text(size = font_size, vjust = 0.35),
          axis.title =   element_text(size= font_size))
  upperPlot
  return(upperPlot)
}


