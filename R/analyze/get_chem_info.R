get_chem_info <- function(all_data, chem_info_old){
  
  all_data <- all_data %>%
    filter(!(chnm %in% c("Tcpp Isomer","Tcpp_isomer")),
           !(chnm == "Chlorpyrifos" & generic_class == "WW"),
           !(chnm == "Caffeine" & generic_class == "WW"),
           !(chnm == "Cotinine" & generic_class == "WW"))
  
  chem_data <-  all_data %>%
    select(SiteID, Date=`Sample Date`, CAS, Value, comment)
  
  chem_info <- select(all_data, CAS, generic_class) %>%
    distinct() %>%
    left_join(select(chem_info_old, CAS, Class, chnm), by="CAS") %>%
    filter(!is.na(CAS)) %>%
    distinct(CAS, .keep_all = TRUE) 
  
  chem_info$Class[is.na(chem_info$Class)] <- chem_info$generic_class[is.na(chem_info$Class)]
  chem_info$Class[chem_info$Class == "pharms"] <- "Pharmaceuticals"
  chem_info <- select(chem_info, -generic_class)
  
  sites_with_detections <- chem_data %>%
    group_by(CAS) %>%
    summarise(n_sites = length(unique(SiteID[Value != 0])))
  
  #DLs:
  dls <- select(all_data, CAS, generic_class, MDL, MQL, Date = `Sample Date`, DL, RL) %>%
    distinct() %>%
    mutate(Date = factor(as.integer(Date), levels = c("2010","2014") )) 
  
  dls$MDL[is.na(dls$MDL)] <- dls$DL[is.na(dls$MDL)]
  dls$MQL[is.na(dls$MQL)] <- dls$RL[is.na(dls$MQL)]
  
  dls$MDL[dls$generic_class == "WW"] <- dls$MDL[dls$generic_class == "WW"]/1000
  dls$MQL[dls$generic_class == "WW"] <- dls$MQL[dls$generic_class == "WW"]/1000
  dls$MDL[dls$generic_class == "pharms"] <- dls$MDL[dls$generic_class == "pharms"]/1000
  dls$MQL[dls$generic_class == "pharms"] <- dls$MQL[dls$generic_class == "pharms"]/1000
  dls$MDL[dls$generic_class %in% c("OC-PCB-PBDE","PAHs")] <- dls$MDL[dls$generic_class  %in% c("OC-PCB-PBDE","PAHs")]/1000000
  dls$MQL[dls$generic_class  %in% c("OC-PCB-PBDE","PAHs")] <- dls$MQL[dls$generic_class  %in% c("OC-PCB-PBDE","PAHs")]/1000000
  
  x <- dls %>%
    select(-DL, -RL) %>%
    gather(variable, value, -Date, -CAS, -generic_class) %>%
    unite(temp, Date, variable) %>%
    spread(temp, value) %>%
    select(-generic_class)
  
  x$`2014_MQL`[x$CAS == "1912-24-9"] <- x$`2014_MQL`[x$CAS == "1912-24-9"][!is.na(x$`2014_MQL`[x$CAS == "1912-24-9"])]
  
  x <- x[-duplicated(x$CAS),]
  
  chem_info_more <- chem_info %>%
    left_join(sites_with_detections, by="CAS") %>%
    left_join(x, by="CAS")
  
  return(chem_info_more)
}

get_exclude <- function(exclude_download){
  exclude <- read.csv(exclude_download, stringsAsFactors = FALSE)
  exclude <- left_join(exclude, select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by = "CAS")
  exclude <- select(exclude, CAS, endPoint, chnm, everything(), -X)
  return(exclude)
}