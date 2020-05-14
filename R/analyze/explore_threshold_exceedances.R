
library(toxEval)
library(tidyverse)
library(readxl)

# Read tox_list and generate chemicalSummary
source("read_chemicalSummary.R")
chem_data <- tox_list$chem_data
chem_site <- tox_list$chem_site
chem_info <- tox_list$chem_info

source("R/analyze/open_land_use.R")
lu <- open_land_use()

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

# Sites monitored each year. Sites in common and sites unique to each year

chem_2010 <- chem_data %>%
  filter(`Sample Date` == 2010) 
sites_2010 <- unique(chem_2010$SiteID)

chem_2014 <- chem_data %>%
  filter(`Sample Date` == 2014) 
sites_2014 <- unique(chem_2014$SiteID)

sum(!(sites_2010 %in% sites_2014))# 32 unique to 2010
sum(!(sites_2014 %in% sites_2010))# 14 unique to 2014
sum((sites_2014 %in% sites_2010))# 23 in common

32+14+23

## How many chemicals ##


#Chems monitored
length(unique(chem_data$CAS))

conc_no_zeros <- chem_data %>% filter(Value > 0)
detected <- unique(conc_no_zeros$CAS)
length(detected)

#number of PAHs detected per site
PAHs <- chem_data %>%
  left_join(chem_info) %>%
  filter(Class == "PAHs", Value > 0) %>%
  group_by(SiteID) %>%
  summarize(num_PAHs = length(unique(CAS)))
mean(PAHs$num_PAHs)
range(PAHs$num_PAHs)

#mean PAH concentrations over all sites
PAH_means <- chem_data %>%
  filter(Value > 0) %>%
  left_join(chem_info) %>%
  filter(Class == "PAHs") %>%
  group_by(CAS,chnm,SiteID) %>%
  summarize(mean_conc = mean(Value)) %>%
  group_by(CAS,chnm) %>%
  summarize(mean_conc = mean(mean_conc)) %>%
  arrange(desc(mean_conc))

# OC pesticides

#mean OC pest concentrations over all sites
OCP_means <- chem_data %>%
  filter(Value > 0) %>%
  left_join(chem_info) %>%
  filter(Class == "OC Pesticides") %>%
  group_by(CAS,chnm,SiteID) %>%
  summarize(mean_conc = mean(Value)) %>%
  group_by(CAS,chnm) %>%
  summarize(mean_conc = mean(mean_conc)) %>%
  arrange(desc(mean_conc))

# Determine the number of chems detected at each site
# Pair this with land use and compute mean number of
# chems for urban, ag, and undeveloped sites

detected_chems_by_site <- chem_data %>%
  left_join(chem_info) %>%
  filter(Value > 0) %>%
  group_by(SiteID) %>%
  summarize(num_chems = length(unique(CAS))) %>%
  left_join(lu, by=c("SiteID"="site"))

mean(detected_chems_by_site$num_chems)
range(detected_chems_by_site$num_chems)
  
urban_thresh <- 15
Ag_thresh <- 40
detected_chems_by_site <- detected_chems_by_site %>%
  mutate(Class_urban = ifelse(Urban > urban_thresh,1,0)) %>%
  mutate(Class_Ag = ifelse(Agriculture > Ag_thresh,1,0)) %>%
  mutate(Class_undev = ifelse(Class_Ag + Class_urban == 0,1,0))

mean_urban_detections <- detected_chems_by_site %>%
  filter(Class_urban == 1) %>%
  summarize(mean_num_chems = mean(num_chems),
            num_sites = length(num_chems))

mean_Ag_detections <- detected_chems_by_site %>%
  filter(Class_Ag == 1) %>%
  summarize(mean_num_chems = mean(num_chems),
            num_sites = length(num_chems))

mean_undev_detections <- detected_chems_by_site %>%
  filter(Class_undev == 1) %>%
  summarize(mean_num_chems = mean(num_chems),
            num_sites = length(num_chems))

 
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
site_exceed$CAS



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

#Look into Menthol endpoints (only one active, relevant endpoint: ACEA_AR_antagonist_80hr)

test <- filter(chemicalSummary,CAS== "89-78-1")
unique(test$endPoint)
range(test$EAR)

test <- filter(test,EAR > 0.001)
unique(test$site)

ACC <- get_ACC("89-78-1")
grep("AR_",ACC$endPoint,)

#plot ACC distribution for the 17 priority chems

plot_ACC_distribution <- function(CAS) {
  library(toxEval)
  library(tidyverse)
  library(ggplot2)
  
  for (i in CAS) {
    ACC_values <- get_ACC(i) %>%
      arrange(ACC)
    
    chem_name <- unique(ACC_values$chnm)
    plot(ACC_values$ACC,ylab="ACC ug/L",main=paste("Distribution of ACC values for",chem_name),pch=20,col="blue")
    
    }
}

#CAS_nums <- c("51218-45-2","75-05-8","124-48-1", "67-66-3", "307-24-4", "1912-24-9", "122-34-9","75-25-2","1763-23-1")

individual_priority_chems <- site_exceed$CAS
mixture_chem_additions <- c("78-40-0","77-93-0", "126-73-8","56-55-3", "205-99-2")

chem_cas <- c(individual_priority_chems,mixture_chem_additions)
  
filename <- "ACC_distributions.pdf"
pdf(filename)
plot_ACC_distribution(chem_cas)
dev.off()
shell.exec(filename)
