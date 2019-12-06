library(toxEval)

explore_endpoints()


library(toxEval)
library(tidyverse)

#NOTE: Add path to path_to_file!!!
path_to_file <- 'data/clean/passive.xlsx' 
tox_list <- create_toxEval(path_to_file)
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC,
                    flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

how_many_eps_per_chem <-group_by(chemical_summary,site,chnm,) %>%
  summarize(num_eps = length(unique(EAR))) %>%
  group_by(chnm) %>%
  summarize(num_eps = max(num_eps)) %>%
  arrange(as.character(chnm))


##
max_by_EP <-group_by(chemical_summary,site,chnm,endPoint) %>%
  summarize(max_EAR_EP = max(EAR))

thresh <- 0.001
site_thresh <- 5

#Number of sites per endpoint per chemical that exceed thresh
site_exceed_by_EP <- group_by(max_by_EP,chnm,endPoint) %>%
  summarize(num_sites_exceeded = sum(max_EAR_EP > thresh)) %>%
  filter(num_sites_exceeded >= site_thresh)

ggplot(data = site_exceed_by_EP,aes(x=chnm,y=num_sites_exceeded)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = site_exceed_by_EP,aes(x=endPoint,y=num_sites_exceeded)) +
  geom_col() +
  facet_wrap(~chnm) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()



#Number of endpoints that exceed thresh for each chemical
exceed_by_EP <- group_by(max_by_EP,site,chnm) %>%
  summarize(num_sites_exceed = sum(max_EAR_EP >= thresh)) %>%
  filter(num_sites_exceed > 0) %>%
  group_by(chnm) %>%
  tally()






DL <- filter(chemical_summary,chnm == "D-Limonene") %>%
  filter(EAR >= thresh)


plot(as.factor(as.character(site_exceed_by_EP$chnm)),site_exceed_by_EP$num_sites_exceeded,las=2)

plot(exceed_by_EP$
