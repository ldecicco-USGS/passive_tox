# Google drive:

# or this:
# https://drive.google.com/drive/folders/10h1mPtyzhkqJx_ikuKx1zoyRD83FZUgF


# To get raw data:
# Easiest...just go to the google drive and download the 3 files:
# https://drive.google.com/drive/folders/10h1mPtyzhkqJx_ikuKx1zoyRD83FZUgF
# "GLRI passive sampler data update 3-17-16.xlsx"
# "GLRI passive sampler pharmaceutical data 8-23-17.xlsx"
# "Analyte Kow and CAS numbers.xlsx"
# library(googledrive)
# dir.create("rawData",showWarnings = FALSE)
# 
# drive_download("GLRI passive sampler data update 3-17-16.xlsx", overwrite = TRUE,
#                path = "rawData/GLRI passive sampler data update 3-17-16.xlsx")
# 
# drive_download(file = "GLRI passive sampler pharmaceutical data 8-23-17.xlsx", overwrite = TRUE,
#                path = "rawData/GLRI passive sampler pharmaceutical data 8-23-17.xlsx")

library(readxl)
library(tidyr)
library(dplyr)
library(openxlsx)

file_2014 <- "rawData/GLRI passive sampler data update 3-17-16.xlsx"
# pharm_file <- "rawData/GLRI passive sampler pharmaceutical data 8-23-17.xlsx"
file_2010 <- "rawData/Copy of Great Lakes passive sampler data update 10-25-13.xlsx"
file_cas <- "rawData/Analyte Kow and CAS numbers.xlsx"
file_pharm2014 <- "rawData/GLRI passive sampler pharmaceutical data UPDATED 8-13-18.xlsx"
file_WW_2014 <- "rawData/GLRI 2014 passive sampler waste indicator.xlsx"

all_cas <- function(file_cas="rawData/Analyte Kow and CAS numbers.xlsx"){

  cas_data <- data.frame()
  
  for(i in c("OC-PCB-PBDE","PAHs","CERC WW","LC8240","LC8069","More")){
    tab_i <- read_excel(file_cas, sheet = i, skip = 3)
    cas_data <- bind_rows(cas_data, tab_i)
  }
  
  cas_data <- cas_data[!is.na(cas_data$Analyte),]
  
  cas_data <- select(cas_data, chnm=Analyte, CAS=`CAS Number`) %>%
    mutate(chnm = tolower(chnm)) %>%
    distinct()
  
  return(cas_data)
}

generic_file_opener <- function(file_name, n_max, sheet, site_sheet,
                                year, skip = 6, skip_site = 3){
  
  data_wide <- read_excel(file_name,
                         sheet = sheet,
                         skip = skip, n_max = n_max)
  
  site_stuff <- read_excel(file_name,
                          sheet = site_sheet,
                          skip = skip_site)
  
  if("CERC Site #" %in% names(site_stuff)){
    site_stuff <- rename(site_stuff, SiteID=`CERC Site #`)
  } else {
    site_stuff <- site_stuff %>%
      rename(SiteID = `CERC ID`,
             STAID = `USGS Station ID`)
  }
  
  site_stuff$SiteID <- gsub("site ","",site_stuff$SiteID, ignore.case = TRUE)
  
  units <- names(data_wide)[-1:-2]
  if(isTRUE(sum(grepl("pg/L", units)) == length(units))){
    convert <- 1000000
  } else if (isTRUE(sum(grepl("ng/L", units)) == length(units))){
    convert <- 1000
  } else if (isTRUE(sum(grepl("ug/L", units)) == length(units))){
    convert <- 1
  } else {
    stop("Check units!")
  }
  
  if(sheet != "est water concentrations"){
    names_wide <- read_excel(file_name,
                             sheet = sheet,
                             skip = 3, n_max = 1) 
    names(data_wide)[4:length(names(names_wide))] <- names(names_wide)[4:length(names(names_wide))]
    names(data_wide)[1] <- "chnm"
    names(data_wide)[2] <- "MDL"
    names(data_wide)[3] <- "MQL"
    data_long <- data_wide %>%
      gather(SiteID, Value, -chnm, -MDL, -MQL) 
  } else {
    names_wide <- read_excel(file_name,
                             sheet = sheet,
                             skip = 4, n_max = 1)
    sites_start <- which(names(names_wide) == "site 1")
    names(data_wide)[sites_start:length(names(names_wide))] <- names(names_wide)[sites_start:length(names(names_wide))]
    names(data_wide)[1] <- "chnm"
    
    if(sites_start == 4){
      names(data_wide)[2] <- "DL"
      names(data_wide)[3] <- "RL"
      data_long <- data_wide %>%
        gather(SiteID, Value, -chnm, -DL,-RL)
      sheet <- "WW"
    } else {
      names(data_wide)[2] <- "Blank"
      data_long <- data_wide %>%
        gather(SiteID, Value, -chnm, -Blank)
      sheet <- "pharms"
    }
  }
  
  data_long$comment <- ""
  data_long$comment[grep("<",data_long$Value)] <- "<"
  data_long$comment[grep("DNQ",data_long$Value)] <- "DNQ"
  data_long$Value <- gsub("DNQ","",data_long$Value)
  data_long$Value <- gsub("<","",data_long$Value)
  data_long$Value <- gsub("a","",data_long$Value)
  data_long$Value <- gsub("b","",data_long$Value)
  data_long$Value <- gsub("c","",data_long$Value)
  data_long$Value <- gsub(" ","",data_long$Value)
  data_long$Value <- as.numeric(data_long$Value) 
  data_long$Value <- data_long$Value/convert
  data_long$generic_class <- sheet
  data_long$`Sample Date` <- year
  data_long$SiteID <- gsub("site ","",data_long$SiteID, ignore.case = TRUE)
  
  data_long <- filter(data_long, 
                      !(is.na(Value) & comment == ""))
  
  data_long <- data_long %>%
    mutate(chnm = tolower(chnm)) %>%
    left_join(all_cas(), by="chnm") %>%
    mutate(chnm = tools::toTitleCase(chnm)) %>%
    left_join(select(site_stuff, SiteID, STAID), by="SiteID") %>%
    mutate(SiteID = dataRetrieval::zeroPad(STAID, 8)) %>%
    select(-STAID) %>%
    distinct() %>%
    filter(!is.na(chnm),
           CAS != "---" | is.na(CAS),
           CAS != "-" | is.na(CAS))

  if(any(is.na(data_long$CAS))){
    message("Some CAS didn't match up")
  }
  
  return(data_long)
}

#####################################################
# OC-PCB-PBDE 2014
OC_2014 <- generic_file_opener(file_2014, 
                               n_max = 45, 
                               sheet = "OC-PCB-PBDE",
                               site_sheet = "site info",
                               year = 2014)

#####################################################
# PAHs 2014:
PAHs_2014 <- generic_file_opener(file_2014, 
                                 n_max = 33, 
                                 sheet = "PAHs",
                                 site_sheet = "site info",
                                 year = 2014)

#####################################################
# Pharm 2014:
pharm_2014 <- generic_file_opener(file_pharm2014, 
                                  n_max = 41, 
                                  sheet = "est water concentrations",
                                  site_sheet = "PrioritySiteInfo",
                                  year = 2014,
                                  skip = 7, skip_site = 2)
pharm_2014$CAS[pharm_2014$chnm == "Buproprion"] <- "34841-39-9"
pharm_2014$CAS[pharm_2014$chnm == "Nadolol"] <- "42200-33-9"

#####################################################
# PAHs 2010:
PAHs_2010 <- generic_file_opener(file_2010,
                                 n_max = 33,
                                 sheet = "PAHs",
                                 site_sheet = "site info",
                                 year = 2010,
                                 skip_site = 2)

#####################################################
# OC-PCB-PBDE 2010
OC_2010 <- generic_file_opener(file_2010, 
                              n_max = 40, 
                              sheet = "OC-PCB-PBDE",
                              site_sheet = "site info",
                              year = 2010,
                              skip_site = 2)
ignore_totals <- c("Total PCBs","Total Pcbs in Mg/l","Total Oc Pesticides")
OC_2010 <- OC_2010[!(OC_2010$chnm %in% ignore_totals),]

#####################################################
# WW 2010
WW_2010 <- generic_file_opener(file_2010, 
                              n_max = 53, 
                              sheet = "WW",
                              site_sheet = "site info",
                              year = 2010,
                              skip_site = 2)
# "Tris(1,3-dichloro-2-propyl)phosphate (T" didn't match up
# looks to be a typo in the data:

WW_2010$chnm[WW_2010$chnm == "Tris(1,3-Dichloro-2-Propyl)Phosphate (t"] <- "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"
WW_2010$CAS[WW_2010$chnm == "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"] <- "13674-87-8"

#####################################################
# WW 2014
WW_2014 <- generic_file_opener(file_WW_2014, 
                               n_max = 46, 
                               sheet = "est water concentrations",
                               site_sheet = "PrioritySiteInfo",
                               year = 2014,
                               skip = 7,
                               skip_site = 2)

#####################################################
# Pharm 2010
pharm_2010 <- generic_file_opener(file_2010, 
                                  n_max = 44, 
                                  sheet = "pharms",
                                  site_sheet = "site info",
                                  year = 2010,
                                  skip_site = 2)

#####################################################
# Sites:
sites_orig_2014 <- read_excel(file_2014,
                     sheet = "site info",
                     skip = 3) %>%
  rename(site_grouping = Lake,
         `Short Name` = `Station shortname`) %>%
  mutate(SiteID = dataRetrieval::zeroPad(STAID, 8)) %>%
  select(SiteID, site_grouping, `Short Name`) %>%
  filter(!is.na(SiteID))

sites_OWC <- data.table::fread("rawData/sites_from_OWC.txt", sep="\t", colClasses = c("SiteID"="character"))

sites_OWC <- select(data.table::setDF(sites_OWC), SiteID, site_grouping, `Short Name`)

sites_orig <- bind_rows(sites_orig_2014, sites_OWC)
sites_orig <- sites_orig[sites_orig$SiteID != "000-----",]

sites_orig_2010 <- read_excel(file_2010,
                              sheet = "site info",
                              skip = 2) %>%
  select(SiteID = `USGS Station ID`) %>%
  filter(!(SiteID %in% sites_orig$SiteID),
         !is.na(SiteID))

sites_orig <- bind_rows(sites_orig, sites_orig_2010)

full_sites <- dataRetrieval::readNWISsite(sites_orig$SiteID)

sites <- sites_orig %>%
  left_join(select(full_sites, SiteID=site_no, Fullname=station_nm, map_nm,
                   dec_lat=dec_lat_va, dec_lon = dec_long_va), by="SiteID")

sites$`Short Name`[is.na(sites$site_grouping)] <- "Pigeon"
sites$site_grouping[is.na(sites$site_grouping)] <- "Lake Superior"

rm(sites_orig, sites_orig_2010, sites_orig_2014, sites_OWC)

##############################################
# Bring it together:
all_data <- bind_rows(pharm_2010,
                      WW_2010,
                      OC_2010,
                      PAHs_2010,
                      pharm_2014,
                      WW_2014,
                      OC_2014,
                      PAHs_2014)

chem_data <- all_data %>%
  select(SiteID, `Sample Date`, CAS, Value, comment)

chem_info_old <- read.csv("data/chem_classes.csv", stringsAsFactors = FALSE)

chem_info <- select(all_data, CAS, generic_class) %>%
  distinct() %>%
  left_join(chem_info_old, by="CAS")

chem_info$Class[is.na(chem_info$Class)] <- chem_info$generic_class[is.na(chem_info$Class)]
chem_info$Class[chem_info$Class == "pharms"] <- "Pharmaceuticals"

chem_info <- chem_info[!duplicated(chem_info$CAS),]
chem_info <- left_join(chem_info, select(tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by = "CAS")

exclude <- read.csv("data/exclude.csv", stringsAsFactors = FALSE)
exclude <- left_join(exclude, select(tox_chemicals, CAS=Substance_CASRN, chnm=Substance_Name), by = "CAS")
exclude <- select(exclude, CAS, endPoint, chnm, everything(), -X)

############################################
# Remove blanks:
blanks <- which(grepl(pattern = "Blank",sites$`Short Name`))
sites <- sites[-blanks,]
chem_data <- chem_data[chem_data$SiteID %in% sites$SiteID,]

# Remove replicates?
resampled <- which(grepl(pattern = "resampled",sites$`Short Name`))
sites <- sites[-resampled,]
chem_data <- chem_data[chem_data$SiteID %in% sites$SiteID,]

# Get rid of censored data:
chem_data$Value[chem_data$comment != ""] <- 0


dir.create("cleanedData",showWarnings = FALSE)

list_of_datasets <- list("Data" = chem_data, 
                         "Chemicals" = chem_info,
                         "Sites" = sites,
                         "Exclude" = exclude)

write.xlsx(list_of_datasets, file = "cleanedData/passive.xlsx", append=TRUE)


