#Process ECOTOX results for chemicals in toxcast

library(tidyverse)


tox_fw <- readRDS("R/Analyze/Out/ECOTOX_combined.Rds")


#Determine stats for each chem
tox_stats <- tox_fw[,-1] %>%
  group_by(chnm, CAS, Chemical.Name, EffectCategory) %>%
  summarize(min_endpoint = min(value),
            max_endpoint = max(value), 
            median_endpoint = median(value),
            fifth_endpoint = quantile(value, probs = 0.05),
            num_endpoints = length(unique(value))) %>%
  ungroup() %>% 
  full_join(chem_CAS) %>%
  select("Class","chnm","CAS","EffectCategory","min_endpoint","fifth_endpoint","median_endpoint", "max_endpoint","num_endpoints","sites_tested","sites_det") %>%
  arrange(is.na(num_endpoints),Class,chnm) %>% 
  rename(`max endpoint` = max_endpoint,
         `5th endpoint percentile` = fifth_endpoint,
         `min endpoint` = min_endpoint)

tox_stats$num_endpoints <- ifelse(is.na(tox_stats$num_endpoints),0,tox_stats$num_endpoints)

tox_stats_ECOTOX <- filter(tox_stats, num_endpoints > 0) %>%
  arrange(EffectCategory,Class,chnm)
max(tox_stats_ECOTOX$num_endpoints)
min(tox_stats_ECOTOX$num_endpoints)
table(tox_stats_ECOTOX$EffectCategory)             # number of chems in ECOTOX + ToxCast


write.csv(tox_stats_ECOTOX,file = "R/Analyze/Out/Tox_endpoint_stats_toxcast.csv",row.names = FALSE)
saveRDS(tox_stats_ECOTOX,file = "R/Analyze/Out/Tox_endpoint_stats_toxcast.rds")


