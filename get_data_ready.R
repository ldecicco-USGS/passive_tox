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
pharm_file <- "rawData/GLRI passive sampler pharmaceutical data 8-23-17.xlsx"
file_2010 <- "rawData/Copy of Great Lakes passive sampler data update 10-25-13.xlsx"
file_cas <- "rawData/Analyte Kow and CAS numbers.xlsx"

generic_file_opener <- function(file_name, n_max, sheet, year, file_cas=file_cas, skip = 6){
  
  data_wide <- read_excel(file_name,
                         sheet = sheet,
                         skip = skip, n_max = n_max)
  
  if(sheet %in% c("pharms","est water concentrations")){
    cas_data_1 <- read_excel(file_cas, sheet = "LC8240", skip = 3)
    cas_data_2 <- read_excel(file_cas, sheet = "LC8069", skip = 3)
    cas_data <- bind_rows(cas_data_1, cas_data_2)
  } else {
    cas_sheet <- switch(sheet,
                        "OC-PCB-PBDE" = "OC-PCB-PBDE",
                        "PAHs" = "PAHs",
                        "WW" = "CERC WW")
    cas_data <- read_excel(file_cas, sheet = cas_sheet, skip = 3)  
  }
  cas_data <- cas_data[!is.na(cas_data$Analyte),]
  
  cas_data <- select(cas_data, chnm=Analyte, CAS=`CAS Number`) %>%
    mutate(chnm = tolower(chnm))
  
  units <- names(data_wide)[-1:-2]
  if(all(grep("pg/L", units) %in% 1:length(units))){
    convert <- 1000000
  } else if (all(grep("ng/L", units) %in% 1:length(units))){
    convert <- 1000
  } else if (all(grep("ug/L", units) %in% 1:length(units))){
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
                             skip = 6, n_max = 1)   
    names(data_wide)[3:length(names(names_wide))] <- names(names_wide)[3:length(names(names_wide))]
    names(data_wide)[1] <- "chnm"
    names(data_wide)[2] <- "Blank"
    data_long <- data_wide %>%
      gather(SiteID, Value, -chnm, -Blank)
    data_long$SiteID <- tools::toTitleCase(data_long$SiteID)
    sheet <- "pharms"
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
  
  data_long <- filter(data_long, 
                      !(is.na(Value) & comment == ""))
  
  data_long <- data_long %>%
    mutate(chnm = tolower(chnm)) %>%
    left_join(cas_data, by="chnm") %>%
    mutate(chnm = tools::toTitleCase(chnm))
  
  if(any(is.na(data_long$CAS))){
    message("Some CAS didn't match up")
  }
  
  
  return(data_long)
}

#####################################################
# OC-PCB-PBDE 2014
data_2014_OC <- generic_file_opener(file_2014, 
                                   n_max = 45, 
                                   sheet = "OC-PCB-PBDE",
                                   year = 2014,
                                   file_cas = file_cas)

#####################################################
# PAHs 2014:
data_2014_PAHs <- generic_file_opener(file_2014, 
                                     n_max = 33, 
                                     sheet = "PAHs",
                                     year = 2014,
                                     file_cas = file_cas)
data_2014_PAHs <- distinct(data_2014_PAHs)
#####################################################
# Pharm 2014:
data_2014_pharm <- generic_file_opener(pharm_file, 
                                      n_max = 38, 
                                      sheet = "est water concentrations",
                                      year = 2014,
                                      file_cas = file_cas,
                                      skip = 10)
# "Nadolol"
# "Buproprion"
data_2014_pharm <- data_2014_pharm[data_2014_pharm$CAS != "-",]
data_2014_pharm$chnm[data_2014_pharm$chnm == "Nadolol"] <- "Nadalol"
data_2014_pharm$CAS[data_2014_pharm$chnm == "Nadalol"] <- "42200-33-9"
data_2014_pharm$chnm[data_2014_pharm$chnm == "Buproprion"] <- "Bupropion"
data_2014_pharm$CAS[data_2014_pharm$chnm == "Bupropion"] <- "34911-55-2"
data_2014_pharm <- data_2014_pharm[!is.na(data_2014_pharm$chnm),]
data_2014_pharm <- distinct(data_2014_pharm)
#####################################################
# PAHs 2010:
data_2010_PAHs <- generic_file_opener(file_2010,
                                      n_max = 33,
                                      sheet = "PAHs",
                                      year = 2010,
                                      file_cas = file_cas)
# "Benzo[g,h,i]perylene" didn't get a CAS
data_2010_PAHs$chnm[data_2010_PAHs$chnm == "Benzo[g,h,i]perylene"] <- "Benzo[g,h,I]perylene"
data_2010_PAHs$CAS[data_2010_PAHs$chnm == "Benzo[g,h,I]perylene"] <- "191-24-2"

#####################################################
# OC-PCB-PBDE 2010
data_2010_OC <- generic_file_opener(file_2010, 
                                    n_max = 40, 
                                    sheet = "OC-PCB-PBDE",
                                    year = 2010,
                                    file_cas = file_cas)
ignore_totals <- c("Total PCBs","Total PCBs in mg/L","Total OC Pesticides")
data_2010_OC <- data_2010_OC[!(data_2010_OC$chnm %in% ignore_totals),]

#####################################################
# WW 2010
data_2010_WW <- generic_file_opener(file_2010, 
                                    n_max = 53, 
                                    sheet = "WW",
                                    year = 2010,
                                    file_cas = file_cas)
# "Tris(1,3-dichloro-2-propyl)phosphate (T" didn't match up
# looks to be a typo in the data:

data_2010_WW$chnm[data_2010_WW$chnm == "Tris(1,3-Dichloro-2-Propyl)Phosphate (t"] <- "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"
data_2010_WW$CAS[data_2010_WW$chnm == "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"] <- "13674-87-8"

data_2010_WW <- data_2010_WW[data_2010_WW$CAS != "---",]
#####################################################
# Pharm 2010
data_2010_Pharm <- generic_file_opener(file_2010, 
                                    n_max = 44, 
                                    sheet = "pharms",
                                    year = 2010,
                                    file_cas = file_cas)

#####################################################
# Sites:
sites_orig <- read_excel(file_2014,
                     sheet = "site info",
                     skip = 3) 

sites_orig <- sites_orig %>%
  rename(SiteID = `CERC Site #`,
         site_grouping = Lake,
         `Short Name` = `Station shortname`) %>%
  mutate(SiteID = paste("Site",SiteID),
         STAID = dataRetrieval::zeroPad(STAID, 8)) %>%
  select(SiteID, site_grouping, `Short Name`, STAID)

full_sites <- dataRetrieval::readNWISsite(sites_orig$STAID[1:49])

sites <- sites_orig %>%
  left_join(select(full_sites, STAID=site_no, Fullname=station_nm, 
                   dec_lat=dec_lat_va, dec_lon = dec_long_va), by="STAID")

site_na <- dataRetrieval::zeroPad(sites$STAID[is.na(sites$dec_lon)],9)
site_na <- site_na[substr(site_na, start = 1, stop = 2) != "00"]
names(site_na) <- substr(site_na, start=2, stop = nchar(site_na))
more_sites <- dataRetrieval::readNWISsite(site_na)

sites_orig$STAID[sites_orig$STAID %in% names(site_na)] <- site_na

sites <- sites_orig %>%
  left_join(select(bind_rows(full_sites, more_sites), 
                   STAID=site_no, Fullname=station_nm, 
                   dec_lat=dec_lat_va, dec_lon = dec_long_va), by="STAID")

sites <- filter(sites, STAID != "000-----")


##############################################
# Bring it together:
all_data <- bind_rows(data_2010_Pharm,
                      data_2010_WW,
                      data_2010_OC,
                      data_2010_PAHs,
                      data_2014_OC,
                      data_2014_PAHs)

chem_data <- all_data %>%
  select(SiteID, `Sample Date`, CAS, Value, comment)

chem_info_old <- read.csv("data/chem_classes.csv", stringsAsFactors = FALSE)

chem_info <- select(all_data, CAS, generic_class) %>%
  distinct() %>%
  left_join(chem_info_old, by="CAS")

chem_info$Class[is.na(chem_info$Class)] <- chem_info$generic_class[is.na(chem_info$Class)]

exclude <- read.csv("data/exclude.csv", stringsAsFactors = FALSE)

############################################
# Remove blanks:
blanks <- which(grepl(pattern = "Blank",sites$`Short Name`))
sites <- sites[-blanks,]
chem_data <- chem_data[chem_data$SiteID %in% sites$SiteID,]

# Remove replicates?
resampled <- which(grepl(pattern = "resampled",sites$`Short Name`))
sites <- sites[-resampled,]
chem_data <- chem_data[chem_data$SiteID %in% sites$SiteID,]


dir.create("cleanedData",showWarnings = FALSE)

list_of_datasets <- list("Data" = chem_data, 
                         "Chemicals" = chem_info,
                         "Sites" = sites,
                         "Exclude" = exclude)

write.xlsx(list_of_datasets, file = "cleanedData/passive.xlsx", append=TRUE)


