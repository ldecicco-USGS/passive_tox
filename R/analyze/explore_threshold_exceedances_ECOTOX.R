
library(toxEval)
library(tidyverse)
library(readxl)

# Read tox_list and generate chemicalSummary

#### Setup ####
path_to_data <- Sys.getenv("PASSIVE_PATH")
# path_to_file_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx")
# path_to_file_non_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx")

path_to_file <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_all.xlsx")
tox_list_ecotox <- create_toxEval(path_to_file)
chemical_summary_ecotox <- get_chemical_summary(tox_list_ecotox)

#Determine number of sites with threshold exceedances of 0.1

## Determine how many chems have EARchem > 0.001 for at least XX % of sites
thresh <- 0.1
Site_proportion_threshold <- 0.1

num_sites_monitored <- chemical_summary_ecotox %>%
  group_by(Class,Bio_category,chnm,CAS) %>%
  summarize(sites_monitored = length(unique(site))) %>%
  arrange(sites_monitored,chnm)

site_exceed_init <- chemical_summary_ecotox %>% 
  group_by(site,Class,chnm,CAS,Bio_category) %>%
  summarize(maxEAR = max(EAR)) %>%
  group_by(Class,chnm, CAS, Bio_category) %>%
  summarize(num_sites_exceeded = sum(maxEAR > thresh)) %>%
  summarize(TQmax = max(maxEAR) %>%
  left_join(num_sites_monitored) %>%
  mutate(proportion_sites_exceeded = num_sites_exceeded/sites_monitored) %>%
  arrange(proportion_sites_exceeded<=Site_proportion_threshold,as.character(chnm))

chem_order <- site_exceed_init %>%
  group_by(CAS,chnm) %>%
  summarize(max_site_exceed = max(proportion_sites_exceeded)) %>%
  arrange(max_site_exceed) %>%
  ungroup() %>%
  mutate(CAS = factor(CAS,levels = CAS))

site_exceed_init$CAS <- factor(site_exceed_init$CAS,levels = chem_order$CAS)

site_exceed_init <- site_exceed_init %>%
  arrange(desc(CAS),Bio_category)

site_exceed <- pivot_wider(site_exceed_init, names_from = Bio_category, values_from = c(num_sites_exceeded,proportion_sites_exceeded)) %>%
  arrange(desc(CAS))

write.csv(site_exceed,file="R/Analyze/Out/ECOTOX_site_threshold_exceedances_all.csv",row.names = FALSE)
saveRDS(site_exceed,file="R/Analyze/Out/ECOTOX_site_threshold_exceedances_all.rds")


