#Process ECOTOX results for chemicals in toxcast

library(tidyverse)
library(ggforce)
library(toxEval)
library(readxl)
library(openxlsx)

chem_CAS <- read.csv(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ToxCast","chem_info_toxcast.csv"),stringsAsFactors = FALSE)
ACC <- get_ACC(chem_CAS$CAS)
file_loc <- file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ECOTOX_priority_chems","explore_priority_exceedances_chemicals.xlsx")
file_loc_csv <- file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ECOTOX_priority_chems","explore_priority_exceedances_chemicals.csv")
chem_CAS_priority <- read_xlsx(path = file_loc,sheet = 1)

files <- list.files(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ToxCast"), pattern="*.txt", all.files=FALSE, full.names=FALSE)

chems <- chem_CAS$chnm
chem_CAS$CAS.Number. <- as.numeric(gsub("-","",chem_CAS$CAS))


for (i in 1:length(files)) {
  tox_temp <- read.delim(file = file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ToxCast",files[i]),sep="|", stringsAsFactors = FALSE)
  names(tox_temp) <- sub("X.","",names(tox_temp))
  tox_temp$Publication.Year <- as.numeric(tox_temp$Publication.Year)
  tox_temp$Conc.1.Mean..Standardized.. <- as.numeric(tox_temp$Conc.1.Mean..Standardized..)
  tox_temp$Conc.Min.1..Standardized.. <- as.numeric(tox_temp$Conc.Min.1..Standardized..)
  tox_temp$Observed.Duration.Mean..Days.. <- as.numeric(tox_temp$Observed.Duration.Mean..Days..)
  #tox_temp$Chemical.Purity.Mean.... <- as.numeric(tox_temp$Chemical.Purity.Mean....)
  tox_temp$Organism.Age.Mean. <- as.numeric(tox_temp$Organism.Age.Mean.)
  tox_temp$Conc.1.Max..Standardized.. <- as.numeric(tox_temp$Conc.1.Max..Standardized..)
  if(i == 1) {tox <- tox_temp 
  }else tox <- bind_rows(tox,tox_temp)
}

#unique(tox$CAS.Number.)
tox <- left_join(tox,chem_CAS)

conc_mean <- ifelse(!(tox$Conc.1.Mean..Standardized.. > 0),NA, tox$Conc.1.Mean..Standardized..)
conc_min <- ifelse(!(tox$Conc.Min.1..Standardized.. > 0),NA, tox$Conc.Min.1..Standardized..)
conc_max <- ifelse(!(tox$Conc.1.Max..Standardized.. > 0),NA, tox$Conc.1.Max..Standardized..)

tox <- transform(tox,value = pmin(conc_mean,conc_min,conc_max,na.rm=TRUE)*1000)
sum(is.na(tox$value))

exposure.type.keep <- c("Aquatic - not reported","Static","Flow-through", "Renewal","Lentic","Lotic")

remove_references <- c(67566, 168095, 171681, 168095, 168095, 11170, 11628)

tox_fw <- tox %>%
  filter(Media.Type == "Fresh water",
         Effect != "Accumulation",
         Exposure.Type %in% exposure.type.keep,
         Conc.1.Type..Standardized..  %in% c("Active ingredient","Total","Formulation"),
         grepl("mg/L",Conc.1.Units..Standardized..),
         !grepl("No significance",Statistical.Significance.),
         !(value < 4.80E-9 & CAS == "1912-24-9"), #Atrazine outlier
         !(value < 0.062 & CAS == "21145-77-7"),  # Tonalide outlier
         !(Reference.Number. %in% remove_references))  #Pyrene study with sediment and pore water


tox_fw_test <- tox_fw %>%  #remove ref 166555 for endpoints that were not adverse effects
  filter((Reference.Number. == 166555 & value < 0.0011))


range(tox_fw_test$value)



# table(tox_fw$Effect)
# unique(tox_fw$Effect.Measurement)
# unique(tox_fw$Exposure.Type)
# table(tox_fw$Conc.1.Units..Standardized..)
# unique(tox_fw$Media.Type)
# table(tox_fw$Statistical.Significance.)
# table(tox_fw$Conc.1.Type..Standardized..)


tox_fw <- tox_fw %>%
  arrange(chnm,value) 

num_chems <- numeric()
chem_index <- numeric()
for (i in 1:length(unique(tox_fw$chnm))) {
  num_chems <- sum(tox_fw$chnm == unique(tox_fw$chnm)[i])
  chem_index <- c(chem_index,1:num_chems)
}
tox_fw$index <- chem_index

benchmark_tab <- tox_fw[,c("CAS.Number.","Chemical.Name","value", "Observed.Duration.Mean..Days..", "Endpoint","Effect","Effect.Measurement")]
names(benchmark_tab) <- c("CAS.Number.","Chemical.Name","Value", "duration", "Endpoint_type","Effect","Effect.Measurement")


#Add Dieldrin benchmark
# dieldrin <- data.frame(60571,"Dieldrin",0.0019,1,"Ambient WQC","","",stringsAsFactors = FALSE)
# names(dieldrin) <- names(benchmark_tab)
# benchmark_tab <- bind_rows(benchmark_tab,dieldrin)

# #Add DDT benchmark
# DDT <- data.frame(50293,"p,p DDT",0.001,1,"Ambient WQC","","",stringsAsFactors = FALSE)
# names(DDT) <- names(benchmark_tab)
# benchmark_tab <- bind_rows(benchmark_tab,DDT)

benchmark_tab <- benchmark_tab %>%
  mutate(endPoint = ifelse(duration > 4,"Chronic","Acute")) %>%
  mutate(groupCol = "ECOTOX")

benchmark_tab <- left_join(benchmark_tab,chem_CAS[,c("CAS.Number.", "CAS","chnm")]) %>%
  rename(Chemical = chnm)

path_to_data <- Sys.getenv("PASSIVE_PATH")

wb <- loadWorkbook(file.path(path_to_data, "data", "data_for_git_repo","clean", "passive.xlsx"))
addWorksheet(wb,sheetName = "Benchmarks")
writeData(wb,sheet = "Benchmarks",x=benchmark_tab)
saveWorkbook(wb,file=file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx"),overwrite = TRUE)
dir.create("R/Analyze/Out", showWarnings = FALSE)
write.csv(tox_fw,"R/Analyze/Out/ECOTOX_filtered_toxcast.csv",row.names = FALSE)
saveRDS(tox_fw,"R/Analyze/Out/ECOTOX_filtered_toxcast.Rds")
# 
# atrazine <- filter(tox_fw,chnm == "Atrazine")
# atrazine$group <- atrazine$Species.Group
# atrazine[grep("Fish",atrazine$group,ignore.case = TRUE),"group"] <- "Fish"
# atrazine[grep("Algae",atrazine$group,ignore.case = TRUE),"group"] <- "Algae"
# atrazine[grep("Insects",atrazine$group,ignore.case = TRUE),"group"] <- "Insects/Spiders"
# atrazine[grep("amphibians",atrazine$group,ignore.case = TRUE),"group"] <- "Amphibians"
# atrazine[grep("Crustaceans",atrazine$group,ignore.case = TRUE),"group"] <- "Crustaceans"
# atrazine[grep("Invertebrates",atrazine$group,ignore.case = TRUE),"group"] <- "Invertebrates"
# atrazine[grep("Molluscs",atrazine$group,ignore.case = TRUE),"group"] <- "Molluscs"
# atrazine[grep("Flowers",atrazine$group,ignore.case = TRUE),"group"] <- "Other plants"
# atrazine[grep("Molluscs",atrazine$group,ignore.case = TRUE),"group"] <- "Molluscs"
# atrazine[grep("Worms",atrazine$group,ignore.case = TRUE),"group"] <- "Worms"
# 
# unique(atrazine$group)
# 
# n_obs <- as.numeric(table(atrazine$group))
# 
# p <- ggplot(atrazine,aes(x=group,y=value)) +
#   geom_boxplot()+
#   scale_y_continuous(trans='log10',limits = c(1e-03,1e+13)) + 
#   theme(axis.text.x = element_text(angle = 90)) +
#   annotate("text",x=1:length(n_obs),y=1e+13,label=n_obs)
# 
# p
# 
# atrazine_low <- filter(atrazine,value < 5)
# n_obs <- as.numeric(table(atrazine_low$group))
# 
# p <- ggplot(atrazine_low,aes(x=group,y=value)) +
#   geom_boxplot()+
#   scale_y_continuous(trans='log10',limits = c(1e-03,10)) + 
#   theme(axis.text.x = element_text(angle = 90)) +
#   annotate("text",x=1:length(n_obs),y=10,label=n_obs) +
#   annotate("text",x=0.6,y=10,label="n =") +
#   ylab("Endpoint Concentration ug/L")
# 
# p
# 
# n_obs <- as.numeric(table(atrazine_low$Effect))
# 
# p <- ggplot(atrazine_low,aes(x=Effect,y=value)) +
#   geom_boxplot()+
#   scale_y_continuous(trans='log10',limits = c(1e-03,10)) + 
#   theme(axis.text.x = element_text(angle = 90)) +
#   annotate("text",x=1:length(n_obs),y=10,label=n_obs) +
#   annotate("text",x=0.6,y=10,label="n =") +
#   ylab("Endpoint Concentration ug/L")
# 
# p
# 
# unique(atrazine_low$Reference.Number.)
# 
# write.csv(atrazine_low,"R/Analyze/Out/Atrazine_less_than_1_endpoints.csv",row.names = FALSE)
# 
# test <- benchmark_tab %>%
#   filter(CAS == "1912-24-9")
# min(test$Value)


p <- ggplot(data = tox_fw,aes(x=Effect,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chnm)

#Determine stats for each chem
tox_stats <- tox_fw[,-1] %>%
  group_by(chnm,CAS,Chemical.Name) %>%
  summarize(min_endpoint = min(value),
            median_endpoint = median(value),
            num_endpoints = length(unique(value))) %>%
  full_join(chem_CAS) %>%
  select("Class","chnm","CAS","min_endpoint","median_endpoint","num_endpoints","sites_tested","sites_det") %>%
  arrange(is.na(num_endpoints),Class,chnm)

tox_stats$num_endpoints <- ifelse(is.na(tox_stats$num_endpoints),0,tox_stats$num_endpoints)

tox_stats_ECOTOX <- filter(tox_stats, num_endpoints > 0)
max(tox_stats_ECOTOX$num_endpoints)
min(tox_stats_ECOTOX$num_endpoints)
dim(tox_stats_ECOTOX)[1]             # number of chems in ECOTOX + ToxCast


write.csv(tox_stats,file = "R/Analyze/Out/Tox_endpoint_stats_toxcast.csv")
saveRDS(tox_stats,file = "R/Analyze/Out/Tox_endpoint_stats_toxcast.rds")

c("Fish","Algae","Amphibians","Crustaceans","Flowers","insects/Spiders","Invertebrates","Molluscs")
