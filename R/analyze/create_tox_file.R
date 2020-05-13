create_tox_object <- function(all_data, chem_info, sites, exclude){
  
  chem_data <- all_data %>%
    filter(!(chnm %in% c("Tcpp Isomer","Tcpp_isomer")),
           !(chnm == "Chlorpyrifos" & generic_class == "WW"),
           !(chnm == "Caffeine" & generic_class == "WW"),
           !(chnm == "Cotinine" & generic_class == "WW")) %>%
    select(SiteID, `Sample Date`, CAS, Value, comment, MDL, MQL, `Date Deployed`, `Date Retrieved`) 
  
  sites_ordered <- sites %>% 
    filter(SiteID %in% chem_data$SiteID)
    
  sites_ordered$map_nm <- substr(gsub("Lake ", "", sites_ordered$site_grouping),1,1)
  
  sites_ordered$map_nm <- paste0(sites_ordered$map_nm, 
                                 c(1:sum(sites_ordered$map_nm == "S"),
                                   1:sum(sites_ordered$map_nm == "M"),
                                   1:sum(sites_ordered$map_nm == "H"),
                                   1:sum(sites_ordered$map_nm == "E"),
                                   1:sum(sites_ordered$map_nm == "O")))
  
  tox_list <- list("Data" = chem_data, 
                   "Chemicals" = chem_info,
                   "Sites" = sites_ordered,
                   "Exclude" = exclude)
  return(tox_list)
}