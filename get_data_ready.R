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

#####################################################
# OC-PCB-PBDE 2014
data_2014_OC <- read_excel(file_2014,
                        sheet = "OC-PCB-PBDE",
                        skip = 6, n_max = 45)
names_OC <- read_excel(file_2014,
                       sheet = "OC-PCB-PBDE",
                       skip = 3, n_max = 1)                       

names(data_2014_OC)[4:length(names(names_OC))] <- names(names_OC)[4:length(names(names_OC))]
names(data_2014_OC)[1] <- "chnm"
names(data_2014_OC)[2] <- "MDL"
names(data_2014_OC)[3] <- "MQL"
rm(names_OC)


data_2014_OC_long <- data_2014_OC %>%
  gather(SiteID, Value, -chnm, -MDL, -MQL) 

data_2014_OC_long$comment <- ""
data_2014_OC_long$comment[grep("<",data_2014_OC_long$Value)] <- "<"
data_2014_OC_long$Value <- gsub("<","",data_2014_OC_long$Value)
data_2014_OC_long$Value <- as.numeric(data_2014_OC_long$Value) 
data_2014_OC_long$`Sample Date` <- 2014
#####################################################
# PAHs
data_2014_PAHs <- read_excel(file_2014,
                           sheet = "PAHs",
                           skip = 6, n_max = 33)
names_PAHs <- read_excel(file_2014,
                       sheet = "PAHs",
                       skip = 3, n_max = 1)                       

names(data_2014_PAHs)[4:length(names(names_PAHs))] <- names(names_PAHs)[4:length(names(names_PAHs))]
names(data_2014_PAHs)[1] <- "chnm"
names(data_2014_PAHs)[2] <- "MDL"
names(data_2014_PAHs)[3] <- "MQL"
rm(names_PAHs)

data_2014_PAHs_long <- data_2014_PAHs %>%
  gather(SiteID, Value, -chnm, -MDL, -MQL) 

data_2014_PAHs_long$comment <- ""
data_2014_PAHs_long$comment[grep("<",data_2014_PAHs_long$Value)] <- "<"
data_2014_PAHs_long$Value <- gsub("<","",data_2014_PAHs_long$Value)
data_2014_PAHs_long$Value <- as.numeric(data_2014_PAHs_long$Value) 
data_2014_PAHs_long$`Sample Date` <- 2014
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

