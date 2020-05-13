library(tidyverse)
library(toxEval)
library(ggforce)

#### Setup ####
library(toxEval)
library(tidyverse)
#NOTE: Add path to path_to_file!!!
path_to_data <- Sys.getenv("PASSIVE_PATH")
path_to_file <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks.xlsx")
path_to_file_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx")
path_to_file_non_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx")
path_to_file_EARs <- file.path(path_to_data, "data", "data_for_git_repo","clean", "passive.xlsx")

tox_list <- create_toxEval(path_to_file)
chemical_summary <- get_chemical_summary(tox_list)
tox_list$benchmarks$acute_chronic <- tox_list$benchmarks$endPoint
tox_list$benchmarks$endPoint <- tox_list$benchmarks$Effect
benchmarks <- tox_list$benchmarks
benchmarks$groupCol <- gsub(benchmarks$groupCol, pattern = " ",replacement = "")

chemical_summary <- get_chemical_summary(tox_list) %>%
  filter(EAR > 0) 
  

### End Setup ###

### Explore Tier 1 vs Tier 2 endpoint concentrations ###

bench_tiers <- benchmarks %>%
  group_by(CAS,groupCol) %>%
  summarize(EP_min = min(ACC_value)) %>%
  pivot_wider(names_from = groupCol,values_from = EP_min) %>%
  left_join(tox_list$chem_info[,c("chnm","CAS","Class")])

tier_scatter <- ggplot(bench_tiers,aes(x=Tier1,y=Tier2)) +
  geom_point() +
  scale_y_continuous(trans='log10',breaks = 10^(seq(-4,7,2))) +
  scale_x_continuous(trans='log10',breaks = 10^(seq(-4,7,2))) +
  geom_abline(intercept = 0, slope = 1, linetype = 1, color = "blue", size = 1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 1, color = "blue", size = 1) + 
  ggtitle(label = "ECOTOX Endpoints for detected passive sampler chemicals",
          subtitle = "Apical/Eco-relevant endpoints (Tier 1) compared to other endpoints (Tier 2)")+ 
  geom_smooth(method='lm', formula= y~x,colour = "orange",se = FALSE)

chems <- tox_list$chem_info[which(tox_list$chem_info$chnm %in% unique(bench_tiers$chnm)),] %>%
  arrange(Class,chnm)



bench_tiers <- bench_tiers %>%
  mutate(ratio = Tier1/Tier2) %>%
  arrange(Class,ratio) %>%
  ungroup() %>%
  mutate(chnm = factor(chnm,levels=chnm))


ratio_tier_comparison <- ggplot(bench_tiers,aes(x=ratio,y=chnm,colour = Class)) +
  geom_point() + 
  scale_x_continuous(trans = 'log10') +
  geom_vline(xintercept = 1,linetype="dashed")

ratio_bp_by_class <- ggplot(bench_tiers,aes(x=ratio,y=Class,colour = Class)) +
  geom_boxplot() + 
  scale_x_continuous(trans = 'log10') +
  geom_vline(xintercept = 1,linetype="dashed")  +
  ggtitle("ECOTOX Endpoint Ratios: min(Tier 1)/min(Tier 2)")

bench_tiers[which(is.na(bench_tiers$Class)),]
unique(tox_list$chem_info$Class)

ear_min_class <- chemical_summary %>%
  group_by(Class) %>%
  summarize(EAR_min = min(EAR)) %>%
  arrange(EAR_min) %>%
  mutate(Class = factor(Class,levels = Class))

# ear_min_chnm <- chemical_summary %>%
#   mutate(Class
#   
# left_join(tox_list$chem_info) %>%
# 
#   
# chemical_summary$chnm <- factor(chemical_summary$chnm, levels = levels(bench_tiers$chnm))

classes <- levels(chemical_summary$Class)

chem_summary_max_EAR <- chemical_summary %>%
  group_by(CAS,Class,site,date,Bio_category) %>%
  summarize(EAR_max_by_sample = max(EAR)) %>%
  group_by(CAS,Class,site,Bio_category) %>%
  summarize(EAR_max = max(EAR_max_by_sample)) %>%
  left_join(tox_list$chem_info[,c("CAS","chnm")])

pdf("R/analyze/plots/Tier1_Tier2_EAR_boxplots_by_class.pdf")
for(i in 1:length(classes)) {
chem_sum_sub <- filter(chem_summary_max_EAR,Class == classes[i])

boxplots_tier1_tier2 <- chem_sum_sub %>%
  ggplot(aes(x=EAR_max,y=chnm,colour = Class)) +
  geom_boxplot() +
  scale_x_continuous(trans = 'log10') + 
  facet_grid(rows = vars(Bio_category),cols = vars(Class)) + 
  geom_vline(xintercept = 0.01,linetype = "dashed")

print(boxplots_tier1_tier2)
             
}
dev.off()

# boxplots_tier1_tier2 <- chemical_summary %>%
#   ggplot(aes(x=EAR,y=chnm,colour = Class)) +
#   geom_boxplot() +
#   scale_x_continuous(trans = 'log10') + 
#   facet_grid_paginate(Bio_category ~ Class,page = 1,ncol = 1,nrow=2,scale = "free_y")
#   # facet_wrap_paginate(.~Bio_category:Class,ncol=2,nrow=2,page=1,scales = "free") +
#   # geom_vline(xintercept = 1,linetype = "dashed")
# 
# boxplots_tier1_tier2

test <- chemical_summary[which(is.na(chemical_summary$chnm)),]



