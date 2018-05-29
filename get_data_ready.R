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

file_2014 <- "rawData/GLRI passive sampler data update 3-17-16.xlsx"
pharm_file <- "rawData/GLRI passive sampler pharmaceutical data 8-23-17.xlsx"
file_2010 <- "rawData/Copy of Great Lakes passive sampler data update 10-25-13.xlsx"

generic_file_opener <- function(file_name, n_max, sheet, year){
  
  data_wide <- read_excel(file_name,
                         sheet = sheet,
                         skip = 6, n_max = n_max)
  names_wide <- read_excel(file_name,
                         sheet = sheet,
                         skip = 3, n_max = 1)                       
  
  names(data_wide)[4:length(names(names_wide))] <- names(names_wide)[4:length(names(names_wide))]
  names(data_wide)[1] <- "chnm"
  names(data_wide)[2] <- "MDL"
  names(data_wide)[3] <- "MQL"

  data_long <- data_wide %>%
    gather(SiteID, Value, -chnm, -MDL, -MQL) 
  
  data_long$comment <- ""
  data_long$comment[grep("<",data_long$Value)] <- "<"
  data_long$Value <- gsub("<","",data_long$Value)
  data_long$Value <- gsub("a","",data_long$Value)
  data_long$Value <- gsub("b","",data_long$Value)
  data_long$Value <- gsub("c","",data_long$Value)
  data_long$Value <- as.numeric(data_long$Value) 
  data_long$`Sample Date` <- year
  
  data_long <- filter(data_long, 
                      !(is.na(Value) & comment == ""))
  
  return(data_long)
}

#####################################################
# OC-PCB-PBDE 2014
data_2014_OC <- generic_file_opener(file_2014, 
                                         n_max = 45, 
                                         sheet = "OC-PCB-PBDE",
                                         year = 2014)

#####################################################
# PAHs 2014:
data_2014_PAHs <- generic_file_opener(file_2014, 
                                         n_max = 33, 
                                         sheet = "PAHs",
                                         year = 2014)

#####################################################
# Sites:
sites <- read_excel(file_2014,
                     sheet = "site info",
                     skip = 3) 

sites <- sites %>%
  rename(site = `CERC Site #`,
         site_grouping = Lake,
         `Short Name` = `Station shortname`) %>%
  mutate(site = paste("Site",site),
         STAID = dataRetrieval::zeroPad(STAID, 8)) %>%
  select(site, site_grouping, `Short Name`, STAID)

full_sites <- dataRetrieval::readNWISsite(sites$STAID[1:49])

sites <- sites %>%
  left_join(select(full_sites, STAID=site_no, station_nm, 
                   dec_lat=dec_lat_va, dec_lon = dec_long_va), by="STAID")


#####################################################
# PAHs 2010:
data_2010_PAHs <- generic_file_opener(file_2010,
                                      n_max = 33,
                                      sheet = "PAHs",
                                      year = 2010)

