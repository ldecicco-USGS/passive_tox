# Script to fill in numbers in the manuscript text

library(tidyverse)
library(toxEval)
source(file = "read_chemicalSummary.R")
source(file = "R/mixtures/mix_script.R")
source(file = "R/analyze/open_land_use.R")

#Land use information
df_lu <- open_land_use()


#Ranges of different land uses
range(df_lu$Developed)
range(df_lu$`Planted/Cultivated`)
range(df_lu$Forest)
range(df_lu$Wetland)

#Get some info on ToxCast assays
# Original assays:
ToxCast_ACC <- ToxCast_ACC

tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))

ToxCast_IN_STUDY <- ToxCast_ACC %>% 
  filter(CAS %in% tox_list$chem_info$CAS)

length(unique(ToxCast_IN_STUDY$endPoint)) #Assays available for chems in this study

length(unique(chemicalSummary$endPoint)) #Assays used in this study after filtering

#Number of assays per chemical
num_assays <- ToxCast_IN_STUDY %>%
  group_by(CAS) %>%
  summarize(num_assays = length(unique(endPoint)))

range(num_assays$num_assays)

num_assays <- chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(CAS) %>%
  summarize(num_assays = length(unique(endPoint)))

range(num_assays$num_assays)# 1-57 assays used per chemical

#CAS numbers for detected chemicals (142 chemicals detected)

x <- tox_list$chem_data %>% 
  filter(Value > 0) %>% 
  select(CAS) %>% 
  distinct() %>% 
  pull()

length(x)

y <- ToxCast_ACC %>% 
  filter(CAS %in% x) %>% 
  select(CAS) %>% 
  distinct()

# Number of detected chemicals in Tox with hits:
nrow(y)
chem_data<- tox_list$chem_data
length(unique(chem_data$CAS[chem_data$Value > 0]))
#THIS DOES INCLUDE PCBs!

# x <- graph_chem_data(chemicalSummary) %>% 
#   filter(meanEAR > 10^-3)
# 
# x %>% 
#   group_by(chnm) %>% 
#   summarize(nsites = length(unique(site))) %>% 
#   filter(nsites >= 10) %>% 
#   ungroup()

ALL_TOX_DATA <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                  "data","data_for_git_repo","raw",
                                  "all_tox_32.rds"))

num_chems_tested <- ALL_TOX_DATA %>% 
  filter(casn %in% x)

# Of those detected, how many are represented in ToxCast (121 chemicals) 
# and how many of those had measureable effects (102)
length(unique(num_chems_tested$casn))
length(unique(num_chems_tested$casn))/length(x)

length(unique(chemicalSummary$chnm[chemicalSummary$EAR > 0]))

chemicalSummary %>% 
!  filter(EAR > 0) %>% 
  group_by(CAS, chnm) %>% 
  summarize(n_eps = length(unique(endPoint))) %>% 
  ungroup() %>% 
  arrange(desc(n_eps))

# How many samples at each site
n_samples <- chemicalSummary %>% 
  select(site, date) %>% 
  distinct() %>%
  group_by(site) %>% 
  summarize(n_samples = length(unique(date))) %>% 
  ungroup() %>% 
  arrange(desc(n_samples)) %>% 
  pull(n_samples) 

#how many sites were sampled (69), and how many of those had more than one sample (24)
length(n_samples)
length(n_samples[n_samples > 1])

#Isolate chemicals not in ToxCast
tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))
CAS_in_study <- tox_list$chem_info$CAS

#Determine chems detected without a corresponding ToxCast assay or no active assays
CAS_detected_in_toxcast <- unique(chemicalSummary$CAS[chemicalSummary$EAR > 0])
CAS_detected <- x
chems_detected_not_in_ToxCast <- x[!(x %in% CAS_detected_in_toxcast)]
chem_info <- tox_list$chem_info
detected_no_ToxCast <- left_join(data.frame(CAS = chems_detected_not_in_ToxCast),tox_list$chem_info[1:3])


#identify chemicals with EAR > threshold of 0.001  
#RESULTS: Num sites exceeded at least once = 29
#         Num sites exceeded at 10% of sites = 12
#         All chems with exceedances were monitored at all 69 sites

thresh <- 0.001
Site_proportion_threshold <- 0

max_EAR_chnm <- chemicalSummary %>% 
  group_by(site,CAS,chnm,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  ungroup() %>% 
  group_by(site,CAS,chnm) %>%
  summarize(EARmax = max(EARsum)) %>%
  group_by(chnm,CAS) %>%
  summarize(Num_sites = length(unique(site)),
            num_sites_exceed = sum(EARmax > thresh),
            EARmax = max(EARmax)) %>%
  ungroup() %>% 
  arrange(desc(EARmax))

exceedances <- max_EAR_chnm %>%
#  filter(EARmax > thresh) %>%
  arrange(desc(num_sites_exceed))

priority_chems_EAR <- exceedances # %>%
#  filter(num_sites_exceed/Num_sites >=0.1)
unique(exceedances$chnm)                                          

saveRDS(priority_chems_EAR,"R/analyze/out/priority_chems_EAR.rds")

num_sites_monitored <- tox_list$chem_data %>% 
  group_by(CAS) %>% 
  summarise(sites_monitored = length(unique(SiteID)))

site_exceed <- chemicalSummary %>% 
  group_by(site,CAS,chnm,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  ungroup() %>% 
  group_by(chnm, CAS, site) %>%
  summarize(maxEAR = max(EARsum),
            num_sites_exceeded = sum(maxEAR > thresh)) %>%
  ungroup() %>% 
  left_join(num_sites_monitored, by = "CAS") %>%
  mutate(proportion_sites_exceeded = num_sites_exceeded/sites_monitored) %>%
  filter(proportion_sites_exceeded > Site_proportion_threshold)



#Menthol comparisons
#Pblication on in-vitro: https://doi.org/10.1016/0009-2797(83)90031-5
#convert mM to ug/L

menthol_endpoints <- c(0.32,0.76) #mM (m-moles/L)
menthol_mw <- 156.27/1000 #g/mole/1000 = g/m-moles
menthol_endpoints_gL <- menthol_endpoints*(menthol_mw)  #g/L
menthol_endpoints_ugL <- menthol_endpoints * 1000000 #ug/L


#ECOTOX text numbers
source(file.path("R","report","chem_priority_summary_Table2.R"))

t2 <- get_table_2()
priority_chems.orig <- readRDS("R/analyze/out/priority_chem_EAR_TQ.rds")

site_thresh <- 0.1

t2$ToxCast <- as.numeric(t2$ToxCast)
t2$ECOTOX_group_1 <- as.numeric(t2$ECOTOX_group_1)
t2$ECOTOX_group_2 <- as.numeric(t2$ECOTOX_group_2)

t2_exceed <- t2 %>%
  mutate(g1_exceed = ECOTOX_group_1/sites_monitored,
         g2_exceed = ECOTOX_group_2/sites_monitored,
         EAR_exceed = ToxCast/sites_monitored,
         max_ecotox = pmax(g1_exceed,g2_exceed,na.rm=TRUE)) %>%
  mutate(g1_boolean = g1_exceed > site_thresh,
         g2_boolean = g2_exceed > site_thresh,
         EAR_boolean = EAR_exceed > site_thresh)

t2_exceed$EAR_G1_exceed <- rowSums(t2_exceed[,c("g1_boolean","EAR_boolean")],na.rm=TRUE)
t2_exceed$EAR_G2_exceed <- rowSums(t2_exceed[,c("EAR_boolean","g2_boolean")],na.rm=TRUE)
t2_exceed$G1_G2_exceed <- rowSums(t2_exceed[,c("g1_boolean","g2_boolean")],na.rm=TRUE)
t2_exceed$sum_exceed <- rowSums(t2_exceed[,c("g1_boolean","g2_boolean","EAR_boolean")],na.rm=TRUE)

#Determine info for text
# 1. how many EAR priorities and which chems
# 2. how many ECOTOX priorities and which chems
# 3. how many Group 1 priorities and which chems
# 4. how many Group 2 priorities and which chems
# 5. What chems match between all three
#                             EAR and Group 1
#                             EAR and Group 2
#                             Group 1 and Group 2

#1. how many EAR priorities and which chems
sum(t2_exceed$EAR_boolean,na.rm=TRUE) # 10 chems
test <- t2_exceed %>%
  filter(EAR_boolean) %>%
  arrange(Class,Chemicals)
test$Chemicals

#2. how many ECOTOX priorities and which chems
sum(t2_exceed$g1_boolean + t2_exceed$g2_boolean > 0,na.rm=TRUE) # 14 chems
test <- t2_exceed %>%
  filter(g1_boolean | g2_boolean) %>%
  arrange(Class,Chemicals)
test$Chemicals

#3. how many Group 1 priorities and which chems
sum(t2_exceed$g1_boolean,na.rm=TRUE) # 14 chems
test <- t2_exceed %>%
  filter(g1_boolean) %>%
  arrange(Class,Chemicals)
group_1 <- test$Chemicals; group_1

#4. how many Group 2 priorities and which chems
sum(t2_exceed$g2_boolean,na.rm=TRUE) # 14 chems
test <- t2_exceed %>%
  filter(g2_boolean) %>%
  arrange(Class,Chemicals)
group_2 <- test$Chemicals; group_2

# How many common among Group 1 and 2
sum(group_2 %in% group_1)

#5. What chems match between all three
#                             EAR and Group 1
#                             EAR and Group 2
#                             Group 1 and Group 2

# EAR and G1 (same as EAR and G2)
test <- t2_exceed %>%
  filter(EAR_boolean & g1_boolean) %>%
  arrange(Class,Chemicals)
test$Chemicals

# EAR and G1 (same as EAR and G2)
test <- t2_exceed %>%
  filter(EAR_boolean & g2_boolean) %>%
  arrange(Class,Chemicals)
test$Chemicals

# G1 and G2
test <- t2_exceed %>%
  filter(g1_boolean & g2_boolean) %>%
  arrange(Class,Chemicals)
test$Chemicals

test <- t2_exceed %>%
  filter(g1_boolean | g2_boolean) %>%
  arrange(Class,Chemicals)
test$Chemicals # 11 total ECOTOX priority chems
               # 2 chems with EAR exceedance, 6 without EAR exceedance, 3 not in toxcast
sum(test$ToxCast < 7,na.rm=TRUE)
sum(test$ToxCast >= 7,na.rm=TRUE)
sum(is.na(test$ToxCast))



priority_chems <- priority_chems.orig %>%
  mutate(g1_exceed = ECOTOX_group_1/sites_monitored,
         g2_exceed = ECOTOX_group_2/sites_monitored,
         EAR_exceed ) %>%
  arrange(desc(g2_exceed)) %>%
  filter(g1_exceed >= site_thresh | g2_exceed >= site_thresh)

group1_chems <- t2 %>%
  filter(!(ECOTOX_group_1=="--")) %>%
  filter(as.numeric(ECOTOX_group_1)/sites_monitored >= site_thresh)
group1_chems$Chemicals  

group2_chems <- t2 %>%
  filter(!(ECOTOX_group_1=="--")) %>%
  filter(as.numeric(ECOTOX_group_2)/sites_monitored >= site_thresh)
group2_chems$Chemicals  

group2_chems$Chemicals   %in% group1_chems$Chemicals  



which(group1_chems$ECOTOX_group_1 == "--")
