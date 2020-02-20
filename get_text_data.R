# Script to fill in numbers in the manuscript text

library(tidyverse)
library(toxEval)
source(file = "read_chemicalSummary.R")
source(file = "R/mixtures/mix_script.R")
source(file = "R/analyze/open_land_use.R")

#Land use information
df_lu <- open_land_use()


#Ranges of different land uses
range(df_lu$Developed)
range(df_lu$`Planted/Cultivated`)
range(df_lu$Forest)
range(df_lu$Wetland)

#Get some info on ToxCast assays
# Original assays:
ToxCast_ACC <- ToxCast_ACC

tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))

ToxCast_IN_STUDY <- ToxCast_ACC %>% 
  filter(CAS %in% tox_list$chem_info$CAS)

length(unique(ToxCast_IN_STUDY$endPoint)) #Assays available for chems in this study

length(unique(chemicalSummary$endPoint)) #Assays used in this study after filtering

#Number of assays per chemical
num_assays <- ToxCast_IN_STUDY %>%
  group_by(CAS) %>%
  summarize(num_assays = length(unique(endPoint)))

range(num_assays$num_assays)

num_assays <- chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(CAS) %>%
  summarize(num_assays = length(unique(endPoint)))

range(num_assays$num_assays)# 1-57 assays used per chemical

#CAS numbers for detected chemicals (142 chemicals detected)

x <- tox_list$chem_data %>% 
  filter(Value > 0) %>% 
  select(CAS) %>% 
  distinct() %>% 
  pull()
length(x)

y <- ToxCast_ACC %>% 
  filter(CAS %in% x) %>% 
  select(CAS) %>% 
  distinct()
# Number of detected chemicals in Tox with hits:
nrow(y)
length(unique(filter(chemicalSummary,EAR > 0)$chnm))

ALL_TOX_DATA <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                  "data","data_for_git_repo","raw",
                                  "all_tox.rds"))

num_chems_tested <- ALL_TOX_DATA %>% 
  filter(casn %in% x)

# Of those detected, how many are represented in ToxCast (121 chemicals) 
# and how many of those had measureable effects (102)
length(unique(num_chems_tested$casn))
length(unique(num_chems_tested$casn))/length(x)
length(unique(chemicalSummary$CAS[chemicalSummary$EAR > 0]))

chemicalSummary %>% 
  filter(EAR > 0) %>% 
  group_by(CAS, chnm) %>% 
  summarize(n_eps = length(unique(endPoint))) %>% 
  ungroup() %>% 
  arrange(desc(n_eps))

# How many samples at each site
n_samples <- chemicalSummary %>% 
  select(site, date) %>% 
  distinct() %>%
  group_by(site) %>% 
  summarize(n_samples = length(unique(date))) %>% 
  ungroup() %>% 
  arrange(desc(n_samples)) %>% 
  pull(n_samples) 

#how many sites were sampled (69), and how many of those had more than one sample (24)
length(n_samples)
length(n_samples[n_samples > 1])

#Isolate chemicals not in ToxCast
tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))
CAS_in_study <- tox_list$chem_info$CAS

#Determine chems detected without a corresponding ToxCast assay or no active assays
CAS_detected_in_toxcast <- unique(chemicalSummary$CAS[chemicalSummary$EAR > 0])
CAS_detected <- x
chems_detected_not_in_ToxCast <- x[!(x %in% CAS_detected_in_toxcast)]
chem_info <- tox_list$chem_info
detected_no_ToxCast <- left_join(data.frame(CAS = chems_detected_not_in_ToxCast),tox_list$chem_info[1:3])


#identify chemicals with EAR > threshold of 0.001  
#RESULTS: Num sites exceeded at least once = 29
#         Num sites exceeded at 10% of sites = 12
#         All chems with exceedances were monitored at all 69 sites

thresh <- 0.001
Site_proportion_threshold <- 0

max_EAR_chnm <- chemicalSummary %>% 
  group_by(site,CAS,chnm,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,CAS,chnm) %>%
  summarize(EARmax = max(EARsum)) %>%
  group_by(chnm,CAS) %>%
  summarize(Num_sites = length(unique(site)),
            num_sites_exceed = sum(EARmax > thresh),
            EARmax = max(EARmax)) %>%
  arrange(desc(EARmax))

exceedances <- max_EAR_chnm %>%
  filter(EARmax > thresh) %>%
  arrange(desc(num_sites_exceed))

unique(exceedances$chnm)                                          


site_exceed <- chemicalSummary %>% group_by(chnm, CAS,site) %>%
  summarize(num_sites_exceeded = sum(maxEAR > thresh)) %>%
  left_join(num_sites_monitored) %>%
  mutate(proportion_sites_exceeded = num_sites_exceeded/sites_monitored) %>%
  filter(proportion_sites_exceeded > Site_proportion_threshold)

