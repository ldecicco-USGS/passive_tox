# Script to fill in numbers in the manuscript text

library(tidyverse)
library(toxEval)
source(file = "read_chemicalSummary.R")
source(file = "R/mixtures/mix_script.R")

#Land use information
df_lu <- open_land_use()

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
         Wetland = `Woody Wetlands` + `Emergent Herbaceous Wetlands`)

#Ranges of different land uses
range(more_lu_cleaned$Developed)
range(more_lu_cleaned$`Planted/Cultivated`)
range(more_lu_cleaned$Forest)
range(more_lu_cleaned$Wetland)

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

y <- ToxCast_ACC %>% 
  filter(CAS %in% x) %>% 
  select(CAS) %>% 
  distinct()

ALL_TOX_DATA <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                  "data","data_for_git_repo","raw",
                                  "all_tox.rds"))

num_chems_tested <- ALL_TOX_DATA %>% 
  filter(casn %in% x)

# Of those detected, how many are represented in ToxCast (121 chemicals) 
# and how many of those had measureable effects (102)
length(unique(num_chems_tested$casn))
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
