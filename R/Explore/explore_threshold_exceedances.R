library(toxEval)

explore_endpoints()


library(toxEval)
library(tidyverse)

#NOTE: Add path to path_to_file!!!
path_to_file <- 'data/clean/passive.xlsx' 
tox_list <- create_toxEval(path_to_file)
chems <- tox_list$chem_info
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC,
                    flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

#How many endpoints per chemical?
how_many_eps_per_chem <-group_by(chemical_summary,site,chnm,) %>%
  summarize(num_eps = length(unique(EAR))) %>%
  group_by(chnm) %>%
  summarize(num_eps = max(num_eps)) %>%
  arrange(as.character(chnm))


## How many sites were chems present and how many of those had EAR > 0.001?

max_by_EP <-group_by(chemical_summary,site,chnm,endPoint) %>%
  summarize(max_EAR_EP = max(EAR))

thresh <- 0.001
site_thresh <- 5
site_occur_thresh <- 10

#Number of sites per endpoint per chemical that exceed thresh
site_exceed_by_EP <- group_by(max_by_EP,chnm,endPoint) %>%
  summarize(num_sites_exceeded = sum(max_EAR_EP > thresh),
            num_sites_present = sum(max_EAR_EP > 0)) %>%
  filter(num_sites_exceeded >= site_thresh & num_sites_present >= site_occur_thresh)

##!!!! work out chnm name matching. Might need to use CAS !!!###
sum(chems$chnm %in% unique(site_exceed_by_EP$chnm))

ggplot(data = site_exceed_by_EP,aes(x=chnm,y=num_sites_exceeded)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot_title <- paste("Chemicals and endpoints present at", site_occur_thresh,
                    "or more sites with EAR > 0.001 at", site_thresh, "or more sites")
ggplot(data = site_exceed_by_EP,aes(x=endPoint,y=num_sites_exceeded)) +
  geom_col() +
  facet_wrap(~chnm) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() +
  ggtitle(label = plot_title)


#Take a look at endopints and what intended target families they represent. 
site_exceed_by_EP$endPoint %in% x$assay_component_name

x_CellCycle <- toxEval::end_point_info %>% filter(intended_target_family == "cell cycle")
site_exceed_by_EP$endPoint %in% x_CellCycle$assay_component_endpoint_name

x_all <- toxEval::end_point_info 
y <- left_join(site_exceed_by_EP,x_all,by=c("endPoint" = "assay_component_endpoint_name"))
y$intended_target_family
unique(y$intended_target_family)

site_exceed_by_EP$endPoint %in% x_all$assay_component_endpoint_name


y <- left_join(site_exceed_by_EP,x,by=c("endPoint" = "assay_component_endpoint_name"))

y$intended_target_family


#Number of endpoints that exceed thresh for each chemical
exceed_by_EP <- group_by(max_by_EP,site,chnm) %>%
  summarize(num_sites_exceed = sum(max_EAR_EP >= thresh)) %>%
  filter(num_sites_exceed > 0) %>%
  group_by(chnm) %>%
  tally()




# How many chemicals with EAR > 0.001 for manuscript text
chems_with_exceedances <- filter(chemical_summary,EAR >= 0.001)
unique(chems_with_exceedances$chnm)  
length(unique(chems_with_exceedances$chnm)  )
