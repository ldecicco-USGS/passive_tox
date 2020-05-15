# Summarize priority chemicals and analyses

path_to_data <- Sys.getenv("PASSIVE_PATH")

EAR_TQ <- readRDS("R/analyze/out/priority_chem_EAR_TQ.rds")
mixtures <- readRDS("R/mixtures/out/mixtures_table.rds")

mixture_chems <- unique(unlist(mixtures$Chemicals))

site_thresh <- 0.1

EAR_TQ$Individual <- ifelse(EAR_TQ$ToxCast/EAR_TQ$sites_monitored > site_thresh,1,0)
EAR_TQ$Group_1 <- ifelse(EAR_TQ$ECOTOX_group_1/EAR_TQ$sites_monitored > site_thresh,1,0)
EAR_TQ$Group_2 <- ifelse(EAR_TQ$ECOTOX_group_2/EAR_TQ$sites_monitored > site_thresh,1,0)

EAR_TQ_priorities <- EAR_TQ %>%
  filter(Individual > 0 | Group_1 > 0 | Group_2 > 0) %>%
  select(c(Class,chnm,Individual, Group_1,Group_2)) %>%
  arrange(Class,chnm)

write.csv(mixture_chems,file=paste0(path_to_data,"/Tables/mixture_chems.csv"))
write.csv(EAR_TQ_priorities,file=paste0(path_to_data,"/Tables/individual_chems.csv"))

