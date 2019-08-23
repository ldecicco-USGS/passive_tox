prep_site_list <- function(Sites){
  
  lakes_ordered <- c("Lake Superior",
                     "Lake Michigan",
                     "Lake Huron",
                     "Lake Erie",
                     "Lake Ontario")
  
  Sites$site_grouping[!(Sites$site_grouping %in% lakes_ordered)] <- "Lake Erie"
  
  Sites$site_grouping <- factor(Sites$site_grouping,
                                levels=lakes_ordered)
  Sites$`Short Name`[Sites$`Short Name` == "Genesee - resampled"] <- "GeneseeDock"
  sites_ordered <- c("StLouis","Nemadji","WhiteWI","Bad",
                     "Montreal","PresqueIsle","Pigeon","Ontonagon",
                     "Sturgeon","Tahquamenon",
                     "Manistique","Escanaba","Ford","Menominee",
                     "Peshtigo","Oconto","Fox","Manitowoc",
                     "Sheboygan #4","Sheboygan #3","Sheboygan #2","Sheboygan #1",
                     "MilwaukeeMouth","IndianaHC #1","IndianaHC #2",
                     "Burns","StJoseph","PawPaw","Kalamazoo",
                     "GrandMI #1","GrandMI #2","GrandMI #3","GrandMI #4",
                     "Muskegon","WhiteMI","PereMarquette","Manistee",
                     "Indian","Cheboygan","ThunderBay","AuSable",
                     "Rifle","Saginaw","BlackMI","Clinton",
                     "Rouge","HuronMI","Raisin",
                     "Maumee","Portage","Sandusky","HuronOH",
                     "Vermilion","BlackOH","Rocky","Cuyahoga",
                     "GrandOH","Ashtabula","Cattaraugus","Buffalo",
                     "Tonawanda","Genesee","GeneseeDock","Oswego","BlackNY",
                     "Oswegatchie","Grass","Raquette","StRegis")
  
  Sites$`Short Name` <- factor(Sites$`Short Name`,
                               levels = sites_ordered)
  
  return(Sites)
  
}


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
    mutate(SiteID = dataRetrieval::zeroPad(SiteID, 8)) %>%
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
  
  sites_ordered <- prep_site_list(sites)

  sites_ordered <- sites_ordered %>%
    arrange(site_grouping, `Short Name`)
  
  return(sites_ordered)
  
}