
library(toxEval)
library(tidyverse)
library(readxl)

# Read tox_list and generate chemicalSummary
source("read_chemicalSummary.R")
chem_data <- tox_list$chem_data
chem_site <- tox_list$chem_site
chem_info <- tox_list$chem_info

chems_hit_calls <- unique(chemicalSummary$CAS)
chem_detected <- chem_data %>%
  filter(Value > 0)
chem_detected_CAS <- unique(chem_detected$CAS)

chem_non_toxcast <- chem_detected_CAS[which(!(chem_detected_CAS %in% chems_hit_calls))]
chem_toxcast <- chem_detected_CAS[which((chem_detected_CAS %in% chems_hit_calls))]

chem_info_non_toxcast <- chem_info %>%
  filter(CAS %in% chem_non_toxcast) %>%
  arrange(Class,chnm)


write.csv(chem_info_non_toxcast,
          file=file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","non-ToxCast","chem_info_non_toxcast.csv"),row.names = FALSE)

chem_info_toxcast <- chem_info %>%
  filter(CAS %in% chem_toxcast) %>%
  arrange(Class,chnm)


write.csv(chem_info_toxcast,
          file=file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ToxCast","chem_info_toxcast.csv"),row.names = FALSE)


# Second method for double checking:

library(tidyverse)
library(toxEval)
source(file = "read_chemicalSummary.R")
tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))
in_out <- tox_list$chem_info %>%
  filter(sites_det > 0) %>%
  select(CAS, chnm, Class) %>%
  left_join(distinct(select(chemicalSummary, CAS)) %>%
              mutate(inTox = TRUE), by="CAS") %>%
  mutate(inTox = ifelse(is.na(inTox), FALSE, TRUE))

not_in <- in_out %>%
  filter(!inTox) %>%
  select(CAS, chnm, Class)

