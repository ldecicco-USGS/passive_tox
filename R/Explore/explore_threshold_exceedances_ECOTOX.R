
library(toxEval)
library(tidyverse)
library(readxl)

# Read tox_list and generate chemicalSummary

#### Setup ####
library(toxEval)
path_to_data <- Sys.getenv("PASSIVE_PATH")
path_to_file <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx")

tox_list <- create_toxEval(path_to_file)
chemical_summary <- get_chemical_summary(tox_list)


source("read_chemicalSummary.R")
chem_data <- tox_list$chem_data
chem_site <- tox_list$chem_site
chem_info <- tox_list$chem_info

source("R/analyze/open_land_use.R")
lu <- open_land_use()

#Determine number of sites with threshold exceedances of 0.1

## Determine how many chems have EARchem > 0.001 for at least XX % of sites
thresh <- 0.01
Site_proportion_threshold <- 0.1

num_sites_monitored <- chemical_summary %>%
  group_by(Class,chnm,CAS) %>%
  summarize(sites_monitored = length(unique(site))) %>%
  arrange(sites_monitored,chnm)

site_exceed <- chemical_summary %>% 
  group_by(site,Class,chnm,CAS) %>%
  summarize(maxEAR = max(EAR)) %>%
  group_by(Class,chnm, CAS) %>%
  summarize(num_sites_exceeded = sum(maxEAR > thresh)) %>%
  left_join(num_sites_monitored) %>%
  mutate(proportion_sites_exceeded = num_sites_exceeded/sites_monitored) %>%
  filter(proportion_sites_exceeded > Site_proportion_threshold) %>%
#  arrange(desc(proportion_sites_exceeded)) %>%
  arrange(as.character(chnm))
write.csv(site_exceed,file="R/Analyze/Out/ECOTOX_site_threshold_exceedances.csv",row.names = FALSE)



