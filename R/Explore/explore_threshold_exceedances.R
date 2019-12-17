
library(toxEval)
library(tidyverse)

#NOTE: Add path to path_to_file!!!
path_to_file <- 'data/clean/passive.xlsx' 
tox_list <- create_toxEval(path_to_file)
chems <- tox_list$chem_info
conc <- tox_list$chem_data
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC,
                    flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)


## EARchem summations ##

EAR_sums <- chemical_summary %>% 
  group_by(site, date,chnm,CAS) %>%
  summarise(sumEAR = sum(EAR)) %>%
  group_by(site,chnm,CAS) %>% 
  summarize(maxEAR = max(sumEAR)) %>% 
  filter(maxEAR > 0)

EAR_exceedances <- filter(EAR_sums, maxEAR > 0.001)
length(unique(EAR_exceedances$site))
length(unique(EAR_exceedances$chnm)) #32 chemicals have EARchem > 0.001 for at least one sample

## Determine how many chems have EARchem > 0.001 for at least XX sites
thresh <- 0.001
site_thresh <- 10


#Number of sites per endpoint per chemical that exceed thresh
site_exceed <- EAR_sums %>% group_by(chnm, CAS,) %>%
  summarize(num_sites_exceeded = sum(maxEAR > thresh)) %>%
  filter(num_sites_exceeded >= site_thresh)
  # 11 sites have EARchem > 0.001 for 10 or more sites
site_exceed$chnm


unique(chemical_summary$site)

# Why menthol??
menthol <- chemical_summary %>% 
  filter(CAS == "89-78-1")
unique(menthol$endPoint)

## End of EARchem summations ##
#############################################


## Individual endpoint exploration ##

#How many endpoints per chemical?
how_many_eps_per_chem <-group_by(chemical_summary,site,chnm,) %>%
  summarize(num_eps = length(unique(EAR))) %>%
  group_by(chnm) %>%
  summarize(num_eps = max(num_eps)) %>%
  arrange(as.character(chnm))


## How many sites were chems present and how many of those had EAR > 0.001?

max_by_EP <-group_by(chemical_summary,site,chnm,CAS,endPoint) %>%
  summarize(max_EAR_EP = max(EAR))

thresh <- 0.001
site_thresh <- 2
site_occur_thresh <- 2

#Number of sites per endpoint per chemical that exceed thresh
site_exceed_by_EP <- group_by(max_by_EP,chnm,CAS,endPoint) %>%
  summarize(num_sites_exceeded = sum(max_EAR_EP > thresh),
            num_sites_present = sum(max_EAR_EP > 0)) %>%
  filter(num_sites_exceeded >= site_thresh & num_sites_present >= site_occur_thresh)

unique(site_exceed_by_EP$chnm)

##!!!! work out chnm name matching. Might need to use CAS !!!###
sum(unique(as.character(site_exceed_by_EP$CAS)) %in% chems$CAS)
priority_individual_chems <- filter(chems,CAS %in% site_exceed_by_EP$CAS)

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
