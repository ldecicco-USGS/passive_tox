create_tox_object <- function(all_data, chem_info, sites, exclude){
  
  chem_data <- all_data %>%
    select(SiteID, `Sample Date`, CAS, Value, comment) %>%
    filter(SiteID %in% sites$SiteID)
    
  tox_list <- list("Data" = chem_data, 
                           "Chemicals" = chem_info,
                           "Sites" = sites,
                           "Exclude" = exclude)
  return(tox_list)
}