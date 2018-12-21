get_chem_info <- function(all_data, chem_info_old){
  
  chem_data <- all_data %>%
    select(SiteID, `Sample Date`, CAS, Value, comment)
  
  chem_info <- select(all_data, CAS, generic_class) %>%
    distinct() %>%
    left_join(chem_info_old, by="CAS")
  
  chem_info$Class[is.na(chem_info$Class)] <- chem_info$generic_class[is.na(chem_info$Class)]
  chem_info$Class[chem_info$Class == "pharms"] <- "Pharmaceuticals"
  
  chem_info <- chem_info[!duplicated(chem_info$CAS),]
  chem_info <- left_join(chem_info, select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by = "CAS")
  
  return(chem_info)
}

get_exclude <- function(exclude_download){
  exclude <- read.csv(exclude_download, stringsAsFactors = FALSE)
  exclude <- left_join(exclude, select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by = "CAS")
  exclude <- select(exclude, CAS, endPoint, chnm, everything(), -X)
  return(exclude)
}