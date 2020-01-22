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

#CAS numbers for detected chemicals (144 chemicals detected)

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

                                        