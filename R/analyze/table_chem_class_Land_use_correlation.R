#Explore land cover relations with different contaminant classes

Chem_Class_correlation_table <- function() {
  library(toxEval)
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(ggpubr)
  
  # -Classes that are likely to have EAR > 0.001
  #   -Endpoints that are important for these classes
  # -Full mixtures
  #   -Endpoints that are important for full mixtures
  #   -Chemicals/classes that play prominent roles in the EARs from mixtures
  
  df_lu <- read_xlsx(path = file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","raw","GLRItox_summary.xlsx"),sheet = 1,skip=1)
  names(df_lu) <- make.names(names(df_lu))
  
  df_ww_lu_2010 <- data.table::fread(file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","raw","WatershedSummary2010_deployment.csv"), 
                                     colClasses = c("USGS_STAID"="character"))
  
  df_ww_lu_2014 <- data.table::fread(file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","raw","WatershedSummary2014_deployment.csv"), 
                                     colClasses = c("USGS_STAID"="character"))
  
  df_ww_lu_2012 <- data.table::fread(file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","raw","WatershedSummary2012_annualaverage.csv"), 
                                     colClasses = c("USGS_STAID"="character"))
  
  
  
  
  chemicalSummary <- readRDS(file = file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","clean","chemical_summary.rds"))
  
  ########################################################################################
  # 1). Determine EARs summation for this analysis using by endpoints within chemical classes 
  # 2). Find max of those EAR summations
  
  # Determine EAR_sum by site/class/endpoint
  class_EP_summary <- chemicalSummary %>%
    group_by(site,date,Class, endPoint) %>%
    summarize(EAR_sum = sum(EAR)) %>%
    group_by(site,Class) 
  
  
  # Determine max EAR_sum (summed by endpoint/class/site)
  class_EP_max <- class_EP_summary %>% 
    summarize(EAR_max = max(EAR_sum),
              EP_max = endPoint[which.max(EAR_sum)])
  
  # Add land use
  lu_columns <- c("STAID", "Urban.......6","Parking.Lot....","Agriculture.......7","Crops....","Water.......14","Wetlands....","Population.Density..people.km2.","Pasture.Hay....")
  names(lu_columns) <- c("site","Urban","Parking_lot","Agriculture","Crops","Water","Wetlands","Population_density","Pasture_Hay")
  
  class_EP_max <- class_EP_max %>% 
    left_join(df_lu[,lu_columns], by=c("site"="STAID"))
  
  names(class_EP_max)[c(1,5:12)] <- names(lu_columns)

  #Determine exceedance of threshold
  class_EP_max$exceed_thresh <- as.integer(ifelse(class_EP_max$EAR_max >= 0.001,1,0))
  
  ### Choose chemical classes that have at least 5% EAR exceedances
  LU <- c("Urban","Crops","Pasture_Hay")
  class_EP_max_tbl <- tbl_df(class_EP_max[,c(LU,"Class","exceed_thresh")])
  
  classes_exceed <- class_EP_max_tbl %>% group_by(Class) %>%
    summarise(exceed_pct = mean(exceed_thresh)) %>%
    filter(exceed_pct > 0.05)
  
  class_EP_max_tbl <- filter(class_EP_max_tbl,(class_EP_max_tbl$Class %in% classes_exceed$Class))
  
  
  options(scipen = 5)
  urban_signif <- class_EP_max_tbl %>%
    group_by(Class) %>%
    do(w = wilcox.test(Urban~exceed_thresh,data=.,paired=FALSE)) %>%
    summarize(Class, urban_p = round(w$p.value,5))
  
  Crop_signif <- class_EP_max_tbl %>%
    group_by(Class) %>%
    do(w = wilcox.test(Crops~exceed_thresh,data=.,paired=FALSE)) %>%
    summarize(Class, crop_p = round(w$p.value,5))
  
  P_H_signif <- class_EP_max_tbl %>%
    group_by(Class) %>%
    do(w = wilcox.test(Pasture_Hay~exceed_thresh,data=.,paired=FALSE)) %>%
    summarize(Class, P_H_p = round(w$p.value,5))
  
  signif <- left_join(urban_signif,Crop_signif) %>%
    left_join(P_H_signif)
  
  signif_best_landuse <- signif %>% pivot_longer(-Class,names_to = "LU",values_to = "p") %>%
    group_by(Class) %>%
    slice(which.min(p)) %>%
    mutate(significant = ifelse(p <= 0.05,1,0))
  
  
  #Summary correlation tables
  LU_signif <- character()
  
  for(i in 1:dim(signif)[1]) {
    sig_cols <- which(signif[i,-1] <= 0.05) + 1
    if(length(sig_cols) > 0) {
      LU_signif[i] <- paste(names(signif)[sig_cols],collapse = "; ")
    }else{
      LU_signif[i] <- ""
    }
    
  }
  
  signif_land_uses <- data.frame(signif$Class,LU_signif)
  
  LU_signif_table <- signif %>%
    mutate(Urban = ifelse(urban_p <= 0.05,"X","")) %>%
    mutate(Crops = ifelse(crop_p <= 0.05,"X","")) %>%
    mutate(Pasture_and_Hay = ifelse(P_H_p <= 0.05,"X","")) %>%
    select(Class,Urban,Crops,Pasture_and_Hay) %>% 
    filter(!(Urban == "" & Crops == "" & Pasture_and_Hay == ""))
  
  return(LU_signif_table)
}
