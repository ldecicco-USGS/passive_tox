#Explore land cover relations with different individual chemicals

Chem_Individual_correlation_table <- function() {
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
  
  df_lu <- read_xlsx(path = file.path("data","raw","GLRItox_summary.xlsx"),sheet = 1,skip=1)
  names(df_lu) <- make.names(names(df_lu))
  
  path_to_file <- 'passive.xlsx' 
  tox_list <- create_toxEval(file.path("data","clean",path_to_file))
  ACC <- get_ACC(tox_list$chem_info$CAS)
  ACC <- remove_flags(ACC = ACC,
                      flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))
  
  cleaned_ep <- clean_endPoint_info(end_point_info)
  filtered_ep <- filter_groups(cleaned_ep, 
                               groupCol = 'intended_target_family',
                               assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                               remove_groups = c('Background Measurement','Undefined'))
  
  chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)
  
  
  ########################################################################################
  # 1). Determine EARs summation for this analysis using by endpoints within chemical classes 
  # 2). Find max of those EAR summations
  
  # Determine EAR_sum by site/chemical/endpoint
  chnm_EP_summary <- chemical_summary %>%
    group_by(site,date,chnm) %>%
    summarize(EAR_sum = sum(EAR)) %>%
    group_by(site,chnm) 
  
  
  # Determine max EAR_sum (summed by endpoint/chnm/site) and the endpoint for that max
  chnm_EP_max <- chnm_EP_summary %>% 
    summarize(EAR_max = max(EAR_sum))
  
  # Add land use
  lu_columns <- c("STAID", "Urban.......6","Parking.Lot....","Agriculture.......7","Crops....","Water.......14","Wetlands....","Population.Density..people.km2.","Pasture.Hay....")
  names(lu_columns) <- c("site","Urban","Parking_lot","Agriculture","Crops","Water","Wetlands","Population_density","Pasture_Hay")
  chnm_EP_max <- left_join(chnm_EP_max,df_lu[,lu_columns],by=c("site"="STAID"))
  names(chnm_EP_max)[c(1,4:11)] <- names(lu_columns)
  
  
  #Determine exceedance of threshold
  chnm_EP_max$exceed_thresh <- as.integer(ifelse(chnm_EP_max$EAR_max >= 0.001,1,0))
  
  ### Choose chemicals that have at least 5% EAR exceedances
  LU <- names(lu_columns)[-1]
  chnm_EP_max_tbl <- tbl_df(chnm_EP_max[,c(LU,"chnm","exceed_thresh")])
  chnm_exceed <- chnm_EP_max_tbl %>% group_by(chnm) %>%
    summarise(exceed_pct = mean(exceed_thresh)) %>%
    filter(exceed_pct > 0.1)
  
  chnm_EP_max_tbl <- filter(chnm_EP_max_tbl,(chnm_EP_max_tbl$chnm %in% chnm_exceed$chnm))
  
  
  options(scipen = 5)
  
  for (i in 1:length(LU)){
    
    LU_signif <- chnm_EP_max_tbl %>%
      mutate(LU_temp = .data[[LU[i]]]) %>%
      group_by(chnm) %>%
      do(w = wilcox.test(LU_temp~exceed_thresh,data=.,paired=FALSE)) %>%
      summarize(chnm, p = round(w$p.value,5))
    names(LU_signif)[2] <- LU[i]
    
    if(i == 1) {signif_accum <- LU_signif 
    }else{
      signif_accum <- left_join(signif_accum,LU_signif)
    }
    
  }
  
  signif_best_landuse <- signif_accum %>% pivot_longer(-chnm,names_to = "LU",values_to = "p") %>%
    group_by(chnm) %>%
    slice(which.min(p)) %>%
    mutate(significant = ifelse(p <= 0.05,1,0))
  
  
  # #Summary correlation tables
  # LU_signif <- character()
  # for(i in 1:dim(signif_accum)[1]) {
  #   sig_cols <- which(signif_accum[i,-1] <= 0.05) + 1
  #   if(length(sig_cols) > 0) {
  #     LU_signif[i] <- paste(names(signif_accum)[sig_cols],collapse = "; ")
  #   }else{
  #     LU_signif[i] <- ""
  #   }
  #   
  # }
  
  LU_signif_table <- signif_accum
  
  for(i in 2:dim(LU_signif_table)[2]){
    LU_signif_table[,i] <- ifelse(signif_accum[,i] <= 0.05,"X","")
  }
  
  return(LU_signif_table)
}