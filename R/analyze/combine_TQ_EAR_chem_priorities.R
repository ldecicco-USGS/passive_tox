
library(toxEval)
library(tidyverse)

source("read_chemicalSummary.R")
chem_data <- tox_list$chem_data
chem_info <- tox_list$chem_info


EAR_priorities <- readRDS("R/analyze/out/priority_chems_EAR.rds")
TQ_all <- readRDS("R/analyze/out/ECOTOX_site_threshold_exceedances_all.rds")
#replace chnm with chem_info version
EAR_priorities <- EAR_priorities[,-grep("chnm",names(EAR_priorities))]
TQ_all <- TQ_all[,-grep("chnm",names(TQ_all))]

TQ_all <- left_join(TQ_all,chem_info[,c("chnm","Class","CAS")])
EAR_priorities <- left_join(EAR_priorities,chem_info[,c("chnm","Class","CAS")])

TQ_all[grep("117-81-7",TQ_all$CAS),"chnm"]
EAR_priorities[grep("117-81-7",EAR_priorities$CAS),"chnm"]
chem_info[grep("117-81-7",chem_info$CAS),"chnm"]

#Process EAR df
EAR_priorities <- EAR_priorities %>%
  rename(sites_monitored = Num_sites, ToxCast = num_sites_exceed)

#Process TQ df
TQ_priorities <- TQ_all %>%
  transform(max_proportion_exceeded = pmax(proportion_sites_exceeded_1,proportion_sites_exceeded_2)) %>%
#  filter(max_proportion_exceeded >= 0.1) %>%
  rename(ECOTOX_group_1 = num_sites_exceeded_1,ECOTOX_group_2 = num_sites_exceeded_2) 
TQ_priorities$sites_monitored

priority_chems <- full_join(EAR_priorities,TQ_priorities) %>%
  select(c(1:5,7,8,12,13)) %>%
  left_join(chem_info[,c("CAS","Class")]) %>%
  mutate(max_exceed = pmax(ToxCast/sites_monitored, 
                           ECOTOX_group_1/sites_monitored, 
                           ECOTOX_group_2/sites_monitored,na.rm = TRUE)) %>%
  mutate(max_sites_exceeded = pmax(ToxCast, 
                           ECOTOX_group_1, 
                           ECOTOX_group_2,na.rm = TRUE)) %>%
  arrange(desc(max_exceed))

write.csv(priority_chems,"R/analyze/out/priority_chem_EAR_TQ.csv",row.names = FALSE)
saveRDS(priority_chems,"R/analyze/out/priority_chem_EAR_TQ.rds")

#priority_chems$CAS %in% chem_info$CAS


