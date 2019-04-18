get_chem_info <- function(all_data, chem_info_old){
  
  chem_data <- all_data %>%
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
  
  chem_info_more <- left_join(chem_info, sites_with_detections, by="CAS")
  
  return(chem_info_more)
}

get_exclude <- function(exclude_download){
  exclude <- read.csv(exclude_download, stringsAsFactors = FALSE)
  exclude <- left_join(exclude, select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by = "CAS")
  exclude <- select(exclude, CAS, endPoint, chnm, everything(), -X)
  return(exclude)
}