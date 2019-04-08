get_sites_ready <- function(file_2014_download, file_2010_download, sites_OWC){
  
  sites_orig_2014 <- readxl::read_excel(file_2014_download,
             sheet = "site info",
             skip = 3) %>%
    rename(site_grouping = Lake,
           `Short Name` = `Station shortname`) %>%
    mutate(SiteID = dataRetrieval::zeroPad(STAID, 8)) %>%
    select(SiteID, site_grouping, `Short Name`) %>%
    filter(!is.na(SiteID))
  
  sites_orig_2014$SiteID[sites_orig_2014$SiteID == "40860038"] <- "040860038"
  sites_orig_2014$SiteID[sites_orig_2014$SiteID == "40851385"] <- "040851385"
  
  sites_orig <- bind_rows(sites_orig_2014, sites_OWC)
  sites_orig <- sites_orig[sites_orig$SiteID != "000-----",]
  
  sites_orig_2010 <- readxl::read_excel(file_2010_download,
                                sheet = "site info",
                                skip = 2) %>%
    select(SiteID = `USGS Station ID`) %>%
    filter(!(SiteID %in% sites_orig$SiteID),
           !is.na(SiteID))

  sites_orig <- bind_rows(sites_orig, sites_orig_2010) %>%
    filter(!duplicated(SiteID))
  
  sites_orig$SiteID[sites_orig$SiteID == "04085790"] <- "04085721"
  
  full_sites <- dataRetrieval::readNWISsite(sites_orig$SiteID)
  
  sites <- sites_orig %>%
    left_join(select(full_sites, SiteID=site_no, Fullname=station_nm, map_nm,
                     dec_lat=dec_lat_va, dec_lon = dec_long_va), by="SiteID")
  
  sites$`Short Name`[is.na(sites$site_grouping)] <- "Pigeon"
  sites$site_grouping[is.na(sites$site_grouping)] <- "Lake Superior"
  
  # blanks <- which(grepl(pattern = "Blank",sites$`Short Name`))
  # sites <- sites[-blanks,]
  
  # resampled <- which(grepl(pattern = "resampled",sites$`Short Name`))
  # sites <- sites[-resampled,]

  return(sites)
  
}