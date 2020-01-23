
library(toxEval)
library(tidyverse)
library(readxl)

# Read tox_list and generate chemicalSummary
source("read_chemicalSummary.R")
chem_data <- tox_list$chem_data
chem_site <- tox_list$chem_site
chem_info <- tox_list$chem_info

# path_to_file <- 'data/clean/passive.xlsx' 
# tox_list <- create_toxEval(path_to_file)
# chems <- tox_list$chem_info
# conc <- tox_list$chem_data
# ACC <- get_ACC(tox_list$chem_info$CAS)
# ACC <- remove_flags(ACC = ACC,
#                     flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))
# 
# cleaned_ep <- clean_endPoint_info(end_point_info)
# filtered_ep <- filter_groups(cleaned_ep, 
#                              groupCol = 'intended_target_family',
#                              assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
#                              remove_groups = c('Background Measurement','Undefined'))
# 
# chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

## How many chemicals ##

# conc_no_zeros <- conc %>% filter(Value > 0)
# detected <- unique(conc_no_zeros$CAS)

# Determine how many sites monitored and how many sites 
# detected per chemical

sites_detected <- chem_info %>%
  mutate(sites_det_fraction = sites_det/sites_tested) %>%
  select(Class,chnm,CAS,sites_det_fraction,sites_tested,sites_det) 

class_order <- sites_detected %>%
  group_by(Class) %>%
  summarize(max_det = max(sites_det_fraction)) %>%
  arrange(desc(max_det))

sites_detected$Class <- factor(sites_detected$Class,levels = class_order$Class)

sites_detected <- sites_detected %>%
  arrange(Class,desc(sites_det_fraction))

ALL_TOX_DATA <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                  "data","data_for_git_repo","raw",
                                  "all_tox.rds"))

sites_detected <- ALL_TOX_DATA %>% 
  filter(casn %in% sites_detected$CAS) %>%
  group_by(casn) %>%
  summarize(active = max(hitc)) %>%
  right_join(sites_detected,by=c("casn"="CAS"))

sites_detected$active[is.na(sites_detected$active)] <- "No ToxCast" #considering all chemicals including non-detects
sites_detected$active[sites_detected$active==0] <- "Inactive"       #considering all chemicals including non-detects
sites_detected$active[sites_detected$active==1] <- "Active"         #considering all chemicals including non-detects

length(unique(sites_detected$casn))
sum(sites_detected$sites_det > 0)

table(filter(sites_detected, sites_det > 0)[,"active"]) #how many of the detected chems were active
table(sites_detected$active)                            #how many of ALL chems were active

site_fraction_threshold <- 0.1
priority_chem_eval <- sites_detected %>%
  group_by(active,chnm,casn,Class) %>%
  summarize(detected_more_than_0.1 = mean(sites_det_fraction > site_fraction_threshold)) %>%
  arrange(active,desc(detected_more_than_0.1),Class)

sum(sites_detected$active != "No ToxCast") #Num chems in ToxCast
sum(sites_detected$active == "No ToxCast") #Num chems not in ToxCast
sum((sites_detected$active != "No ToxCast") & (sites_detected$sites_det_fraction>0)) #Num chems detected in ToxCast
sum(((sites_detected$active == "No ToxCast") & (sites_detected$sites_det_fraction>0))) #Num chems detected in ToxCast


sum(sites_detected$active == "No ToxCast")
sum(sites_detected$active == "No ToxCast" & sites_detected$sites_det_fraction> site_fraction_threshold)
sites_detected_no_ToxCast <- filter(sites_detected, active == "No ToxCast")

#How many of each class are not in ToxCast but have 10% or more site detections
sites_detected_no_ToxCast %>%
  filter(sites_det_fraction > site_fraction_threshold) %>%
  group_by(Class) %>%
  summarize(num_chems_per_class = length(chnm))
# Of those detected, how many are represented in ToxCast (121 chemicals) 
# and how many of those had measureable effects (102)
length(unique(num_chems_tested$casn))
length(unique(chemicalSummary$CAS[chemicalSummary$EAR > 0]))

#list of chemicals with detections, represented in ToxCast with EAR > 0
casn_in_toxcast_detected <- unique(chemicalSummary$CAS[chemicalSummary$EAR > 0])
casn_in_toxcast_detected2 <- sites_detected[sites_detected$active == "Active" & sites_detected$sites_det>0,]$casn
casn_in_toxcast_detected %in% casn_in_toxcast_detected2
sum(!(casn_in_toxcast_detected2 %in% casn_in_toxcast_detected))

detection_summary <- sites_detected %>%
  mutate(site_thresh_exceed = sites_det_fraction > site_fraction_threshold) %>%
  group_by(active,Class,site_thresh_exceed) %>%
  filter(site_thresh_exceed) %>%
  arrange(active) 

# ADD NUM CHEMS IN TOXCAST AND ACTIVE
# EXPLORE INACTIVE AND NOT IT TOXCAST FOR DETECTION FREQUENCY


## EARchem summations ##

EAR_sums <- chemicalSummary %>% 
  group_by(site, date,chnm,CAS) %>%
  summarise(sumEAR = sum(EAR)) %>%
  group_by(site,chnm,CAS) %>% 
  summarize(maxEAR = max(sumEAR))

num_sites_monitored <- EAR_sums %>%
  group_by(chnm,CAS) %>%
  summarize(sites_monitored = length(unique(site))) %>%
  arrange(sites_monitored,chnm)

EAR_sums <- EAR_sums %>% 
  filter(maxEAR > 0)

EAR_exceedances <- filter(EAR_sums, maxEAR > 0.001)
length(unique(EAR_exceedances$site))
length(unique(EAR_exceedances$chnm)) #32 chemicals have EARchem > 0.001 for at least one sample

## Determine how many chems have EARchem > 0.001 for at least XX % of sites
thresh <- 0.001
Site_proportion_threshold <- 0.1

site_exceed <- EAR_sums %>% group_by(chnm, CAS) %>%
  summarize(num_sites_exceeded = sum(maxEAR > thresh)) %>%
  left_join(num_sites_monitored) %>%
  mutate(proportion_sites_exceeded = num_sites_exceeded/sites_monitored) %>%
  filter(proportion_sites_exceeded > Site_proportion_threshold)

# 12 chemicals have EARchem > 0.001 for 10% or more sites monitored
site_exceed$chnm




## End of priority chemical determination ##
#############################################




## Individual endpoint exploration ##

#How many endpoints per chemical?
how_many_eps_per_chem <-group_by(chemicalSummary,site,chnm,) %>%
  summarize(num_eps = length(unique(EAR))) %>%
  group_by(chnm) %>%
  summarize(num_eps = max(num_eps)) %>%
  arrange(as.character(chnm))


## How many sites were chems present and how many of those had EAR > 0.001?

max_by_EP <-group_by(chemicalSummary,site,chnm,CAS,endPoint) %>%
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
chems_with_exceedances <- filter(chemicalSummary,EAR >= 0.001)
unique(chems_with_exceedances$chnm)  
length(unique(chems_with_exceedances$chnm)  )


#Look into cis vs trans isomers
priority_chem_eval[ grep("cis",priority_chem_eval$chnm,ignore.case = TRUE),]
priority_chem_eval[ grep("trans",priority_chem_eval$chnm,ignore.case = TRUE),]

