
library(toxEval)
library(tidyverse)
library(readxl)

# Read tox_list and generate chemicalSummary

#### Setup ####
library(toxEval)
path_to_data <- Sys.getenv("PASSIVE_PATH")
path_to_file_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx")
path_to_file_non_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx")


tox_list_toxcast <- create_toxEval(path_to_file_toxcast)
chemical_summary_toxcast <- get_chemical_summary(tox_list_toxcast)

tox_list_non_toxcast <- create_toxEval(path_to_file_non_toxcast)
chemical_summary_non_toxcast <- get_chemical_summary(tox_list_non_toxcast)


# chem_data <- tox_list$chem_data
# chem_site <- tox_list$chem_site
# chem_info <- tox_list$chem_info

# source("R/analyze/open_land_use.R")
# lu <- open_land_use()

#Determine number of sites with threshold exceedances of 0.1

## Determine how many chems have EARchem > 0.001 for at least XX % of sites
thresh <- 0.1
Site_proportion_threshold <- 0.1

num_sites_monitored <- chemical_summary_toxcast %>%
  group_by(Class,chnm,CAS) %>%
  summarize(sites_monitored = length(unique(site))) %>%
  arrange(sites_monitored,chnm)

site_exceed_toxcast <- chemical_summary_toxcast %>% 
  group_by(site,Class,chnm,CAS) %>%
  summarize(maxEAR = max(EAR)) %>%
  group_by(Class,chnm, CAS) %>%
  summarize(num_sites_exceeded = sum(maxEAR > thresh)) %>%
  left_join(num_sites_monitored) %>%
  mutate(proportion_sites_exceeded = num_sites_exceeded/sites_monitored) %>%
#  filter(proportion_sites_exceeded > Site_proportion_threshold) %>%
#  arrange(desc(proportion_sites_exceeded))
  arrange(as.character(chnm))
write.csv(site_exceed_toxcast,file="R/Analyze/Out/ECOTOX_site_threshold_exceedances_toxcast.csv",row.names = FALSE)


num_sites_monitored <- chemical_summary_non_toxcast %>%
  group_by(Class,chnm,CAS) %>%
  summarize(sites_monitored = length(unique(site))) %>%
  arrange(sites_monitored,chnm)

site_exceed_non_toxcast <- chemical_summary_non_toxcast %>% 
  group_by(site,Class,chnm,CAS) %>%
  summarize(maxEAR = max(EAR)) %>%
  group_by(Class,chnm, CAS) %>%
  summarize(num_sites_exceeded = sum(maxEAR > thresh)) %>%
  left_join(num_sites_monitored) %>%
  mutate(proportion_sites_exceeded = num_sites_exceeded/sites_monitored) %>%
  #  filter(proportion_sites_exceeded > Site_proportion_threshold) %>%
  arrange(desc(proportion_sites_exceeded))
#  arrange(as.character(chnm))
write.csv(site_exceed,file="R/Analyze/Out/ECOTOX_site_threshold_exceedances_toxcast.csv",row.names = FALSE)



#Test why TCEP does not show up on the priority list

test <- chemical_summary_toxcast %>%
  filter(grepl("TDCPP",chnm),EAR > 0)
boxplot(test$EAR,log="y")
abline(h=0.1)

test <- chemical_summary_toxcast %>%
  group_by(site,date,chnm)%>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,date,chnm) %>%
  summarize(EARsummax = max(EARsum)) %>%
  filter(grepl("TCEP",chnm),EARsummax > 0)
boxplot(test$EARsummax,log="y",horizontal = TRUE)
abline(v=0.1)

test2 <- test %>%
  filter(EARsummax > 0.1)
unique(test2$site)


test <- chemical_summary_toxcast %>%
  group_by(site,date,chnm)%>%
  summarize(EARmax = max(EAR)) %>%
  filter(grepl("TCEP",chnm),EARmax > 0)
boxplot(test$EARmax,log="y",horizontal = TRUE)
abline(v=0.1)
