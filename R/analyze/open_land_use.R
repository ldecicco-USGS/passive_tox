open_land_use <- function(){
  
  not_all_na <- function(x) all(!is.na(x))
  
  df_lu <- readxl::read_xlsx(path = file.path(Sys.getenv("PASSIVE_PATH"),
                                              "data","data_for_git_repo","raw",
                                              "GLRItox_summary.xlsx"),
                             sheet = 1, skip=1) %>% 
    rename(site = STAID, 
           Basin_Area_mi2 = `Basin Area (mi2)`,
           Basin_area_km2 = `Basin Area (km2)`,
           Urban = `Urban (%)...9`, 
           Parking_lot = `Parking Lot (%)`,
           Agriculture = `Agriculture (%)...10`) %>% 
    # mutate(Developed = Urban + Agriculture) %>% 
    select(-`Urban (%)...6`,   # this is getting rid of 2011 data
           -`Agriculture (%)...7`, 
           -`Other (%)...8`,
           -`Impervious Area (%)...12`,
           -`[Outside of dataset extent]`) %>% 
    select_if(not_all_na)
  
  names(df_lu) <- gsub("\\s*\\([^\\)]+\\)",
                       replacement = "",
                       names(df_lu))
  names(df_lu) <- gsub(pattern = "\\...",
                       replacement = "",
                       names(df_lu))
  names(df_lu) <- gsub(pattern = ", ",
                       replacement = "_",
                       names(df_lu))
  names(df_lu) <- gsub(pattern = "/",
                       replacement = "_",
                       names(df_lu))
  names(df_lu) <- gsub(pattern = " ",
                       replacement = "_",
                       names(df_lu))
  names(df_lu) <- gsub(pattern = "\\[",
                       replacement = "",
                       names(df_lu))
  names(df_lu) <- gsub(pattern = "]",
                       replacement = "",
                       names(df_lu))
  names(df_lu) <- gsub(pattern = "-",
                       replacement = "_",
                       names(df_lu))
  
  # open other stuff:
  watershed2010 <- data.table::fread(file.path(Sys.getenv("PASSIVE_PATH"),
                                               "data","data_for_git_repo","raw",
                                               "WatershedSummary2010.csv"), 
                                     colClasses = c("USGS_STAID"="character"),
                                     data.table = FALSE) %>% 
    select(DWfrac_2010 = DW_FlowFraction,
           frac_2010 = FlowFraction,
           effluent_2010 = Effluent_mgalperday,
           site = USGS_STAID)
  
  watershed2014 <- data.table::fread(file.path(Sys.getenv("PASSIVE_PATH"),
                                               "data","data_for_git_repo","raw",
                                               "WatershedSummary2014.csv"), 
                                     colClasses = c("USGS_STAID"="character"),
                                     data.table = FALSE) %>% 
    select(DWfrac_2014 = DW_FlowFraction,
           frac_2014 = FlowFraction,
           effluent_2014 = Effluent_mgalperday,
           site = USGS_STAID)
  
  combo <- watershed2010 %>% 
    left_join(watershed2014, by = "site") 
  
  df_lu_new <- df_lu %>% 
    left_join(combo, by = "site") %>% 
    mutate_if(is.numeric,
              list(~ case_when(is.na(.) ~ 0, 
                               TRUE ~ .) )) %>% 
    filter(site != "04010500")
  
  # Need to normalize the fraction stuff?
  
  df_lu_new$DWfrac_2010 <- 100*df_lu_new$DWfrac_2010
  df_lu_new$DWfrac_2014 <- 100*df_lu_new$DWfrac_2014
  
  more_lu <- readxl::read_xlsx(path = file.path(Sys.getenv("PASSIVE_PATH"),
                                                "data","data_for_git_repo","raw",
                                                "GLRItox_summary.xlsx"),
                               sheet = "NLCD_LC2016", n_max = 184)
  
  more_lu_cleaned <- more_lu %>% 
    select(site =  AREAID,
           `Open Water` = PCT_11,
           `Developed, Open Space` = PCT_21,
           `Developed, Low Intensity` = PCT_22,
           `Developed, Medium Intensity` = PCT_23,
           `Developed High Intensity` = PCT_24,
           `Barron` = PCT_31,
           `Deciduous Forest` = PCT_41,
           `Evergreen Forest` = PCT_42,
           `Mixed Forest` = PCT_43,
           `Shrubland` = PCT_52,
           `Herbaceous` = PCT_71,
           `Pasture/Hay` = PCT_81,
           `Cultivated Crops` = PCT_82,
           `Woody Wetlands` = PCT_90,
           `Emergent Herbaceous Wetlands` = PCT_95) %>% 
    mutate(Developed = `Developed, Open Space` + `Developed, Low Intensity` +
             `Developed, Medium Intensity` + `Developed High Intensity`,
           Forest = `Deciduous Forest` + `Mixed Forest`,
           `Planted/Cultivated` = `Pasture/Hay` + `Cultivated Crops`,
           Wetland = `Woody Wetlands` + `Emergent Herbaceous Wetlands`) %>% 
    select(site, Developed, Forest, Wetland, `Planted/Cultivated`)
  
  df_lu_new <- df_lu_new %>% 
    left_join(more_lu_cleaned, by="site")
  
  return(df_lu_new)
  
}
