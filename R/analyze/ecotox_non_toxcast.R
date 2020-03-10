library(tidyverse)
library(toxEval)
library(readxl)
library(openxlsx)

chem_CAS <- read.csv(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","non-ToxCast","chem_info_non_toxcast.csv"),stringsAsFactors = FALSE)
ACC <- get_ACC(chem_CAS$CAS)

files <- "ECOTOX_nontoxcast.txt"
  #list.files(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","non-ToxCast"), pattern="*.txt", all.files=FALSE, full.names=FALSE)

chems <- chem_CAS$chnm
chem_CAS$CAS.Number. <- as.numeric(gsub("-","",chem_CAS$CAS))


for (i in 1:length(files)) {
  tox_temp <- read.delim(file = file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","non-ToxCast",files[i]),sep="|", stringsAsFactors = FALSE)
  names(tox_temp) <- sub("X.","",names(tox_temp))
  tox_temp$Conc.1.Mean..Standardized.. <- as.numeric(tox_temp$Conc.1.Mean..Standardized..)
  tox_temp$Conc.Min.1..Standardized.. <- as.numeric(tox_temp$Conc.Min.1..Standardized..)
  tox_temp$Observed.Duration.Mean..Days.. <- as.numeric(tox_temp$Observed.Duration.Mean..Days..)
  #tox_temp$Chemical.Purity.Mean.... <- as.numeric(tox_temp$Chemical.Purity.Mean....)
  tox_temp$Organism.Age.Mean. <- as.numeric(tox_temp$Organism.Age.Mean.)
  tox_temp$Conc.1.Max..Standardized.. <- as.numeric(tox_temp$Conc.1.Max..Standardized..)
  if(i == 1) {tox <- tox_temp 
  }else tox <- bind_rows(tox,tox_temp)
}

tox <- left_join(tox,chem_CAS)

conc_mean <- ifelse(!(tox$Conc.1.Mean..Standardized.. > 0),NA, tox$Conc.1.Mean..Standardized..)
conc_min <- ifelse(!(tox$Conc.Min.1..Standardized.. > 0),NA, tox$Conc.Min.1..Standardized..)
conc_max <- ifelse(!(tox$Conc.1.Max..Standardized.. > 0),NA, tox$Conc.1.Max..Standardized..)

tox <- transform(tox,value = pmin(conc_mean,conc_min,conc_max,na.rm=TRUE)*1000)
sum(is.na(tox$value))

exposure.type.keep <- c("Aquatic - not reported","Static","Flow-through", "Renewal","Lentic","Lotic")

tox_fw <- tox %>%
  filter(Media.Type == "Fresh water",
         Effect != "Accumulation",
         Exposure.Type %in% exposure.type.keep,
         Conc.1.Type..Standardized.. == "Active ingredient",
         grepl("mg/L",Conc.1.Units..Standardized..),
         !grepl("No significance",Statistical.Significance.))



# 
# table(tox_fw$Effect)
# unique(tox_fw$Effect.Measurement)
# unique(tox_fw$Exposure.Type)
# table(tox_fw$Conc.1.Units..Standardized..)
# unique(tox_fw$Media.Type)
# table(tox_fw$Statistical.Significance.)
# table(tox_fw$Conc.1.Type..Standardized..)
# 
 

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


#Add PCB benchmark
pcbs <- data.frame(1336363,"Total PCBs",0.0015,1,"Ambient WQC","","",stringsAsFactors = FALSE)
names(pcbs) <- names(benchmark_tab)

benchmark_tab <- bind_rows(benchmark_tab,pcbs)

benchmark_tab <- benchmark_tab %>%
  mutate(endPoint = ifelse(duration > 4,"Chronic","Acute")) %>%
  mutate(groupCol = "ECOTOX")

benchmark_tab <- left_join(benchmark_tab,chem_CAS[,c("CAS.Number.", "CAS","chnm")]) %>%
  rename(Chemical = chnm)
#benchmark_tab <- full_join(benchmark_tab,pest_benchmarks)


path_to_data <- Sys.getenv("PASSIVE_PATH")

wb <- loadWorkbook(file.path(path_to_data, "data", "data_for_git_repo","clean", "passive.xlsx"))
addWorksheet(wb,sheetName = "Benchmarks")
writeData(wb,sheet = "Benchmarks",x=benchmark_tab)
saveWorkbook(wb,file=file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx"),overwrite = TRUE)
saveRDS(tox_fw,"R/Analyze/Out/ECOTOX_filtered_non_toxcast.rds")
write.csv(tox_fw,"R/Analyze/Out/ECOTOX_filtered_non_toxcast.csv",row.names = FALSE)

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

write.csv(tox_stats,file = "R/Analyze/Out/Tox_endpoint_stats_non_toxcast.csv")
saveRDS(tox_stats,file = "R/Analyze/Out/Tox_endpoint_stats_non_toxcast.rds")

#c("Fish","Algae","Amphibians","Crustaceans","Flowers","insects/Spiders","Invertebrates","Molluscs")
