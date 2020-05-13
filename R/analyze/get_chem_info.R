get_chem_info <- function(all_data, chem_info_old){
  
  all_data <- all_data %>%
    filter(!(chnm %in% c("Tcpp Isomer","Tcpp_isomer")),
           !(chnm == "Chlorpyrifos" & generic_class == "WW"),
           !(chnm == "Caffeine" & generic_class == "WW"),
           !(chnm == "Cotinine" & generic_class == "WW"))
  
  chem_data <-  all_data %>%
    select(SiteID, Date=`Sample Date`, CAS, Value, comment)
  
  chem_info <- select(all_data, CAS, generic_class, chnm) %>%
    distinct() %>%
    left_join(distinct(select(chem_info_old, CAS, Class)), by="CAS") %>%
    filter(!is.na(CAS)) %>%
    distinct(CAS, .keep_all = TRUE) 
  
  chem_info$Class[is.na(chem_info$Class)] <- chem_info$generic_class[is.na(chem_info$Class)]
  chem_info$Class[chem_info$Class == "pharms"] <- "Pharmaceuticals"
  chem_info <- select(chem_info, -generic_class)
  
  sites_with_detections <- chem_data %>%
    group_by(CAS) %>%
    summarise(sites_tested = length(unique(SiteID)),
              sites_det = length(unique(SiteID[Value != 0])))
  
  #DLs:
  dls <- select(all_data, CAS, generic_class, MDL, MQL, Date = `Sample Date`) %>% #, DL, RL) %>%
    distinct() 

  x <- dls %>%
    tidyr::gather(variable, value, -Date, -CAS, -generic_class) %>%
    tidyr::unite(temp, Date, variable) %>%
    tidyr::spread(temp, value) %>%
    dplyr::select(-generic_class)

  x$`2014_MQL`[x$CAS == "1912-24-9"] <- x$`2014_MQL`[x$CAS == "1912-24-9"][!is.na(x$`2014_MQL`[x$CAS == "1912-24-9"])]
  
  x <- x[-duplicated(x$CAS),]
  
  chem_info_more <- chem_info %>%
    left_join(sites_with_detections, by="CAS") %>%
    left_join(x, by="CAS")
  
  chem_info_more <- chem_info_more[!duplicated(chem_info_more$CAS, fromLast = TRUE),]
  
  return(chem_info_more)
}

get_exclude <- function(exclude_download){
  exclude <- data.table::fread(exclude_download, data.table = FALSE)
  exclude <- left_join(exclude, 
                       select(toxEval::tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), 
                       by = c("CAS"))
  exclude <- select(exclude, CAS, endPoint, chnm, dplyr::everything())
  return(exclude)
}


fix_cas <- function(df, cas_change){
  
  cas_change_clean <- cas_change %>%
    select(CAS = `Original CAS`, new_CAS = `CAS to update in analytical data`) %>%
    filter(!is.na(CAS))
  
  df_fixed <- df %>%
    left_join(cas_change_clean, by="CAS")
  
  df_fixed$CAS[!is.na(df_fixed$new_CAS)] <- df_fixed$new_CAS[!is.na(df_fixed$new_CAS)]
  
  df_fixed <- select(df_fixed, -new_CAS)

  return(df_fixed)  
}