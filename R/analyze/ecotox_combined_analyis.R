#Process ECOTOX results for chemicals in toxcast

library(tidyverse)
library(ggforce)
library(toxEval)
library(readxl)
library(openxlsx)

chem_CAS_toxcast <- read.csv(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ToxCast","chem_info_toxcast.csv"),stringsAsFactors = FALSE)
chem_CAS_toxcast$toxcast <- "Yes"
chem_CAS_non_toxcast <- read.csv(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","non-ToxCast","chem_info_non_toxcast.csv"),stringsAsFactors = FALSE)
chem_CAS_non_toxcast$toxcast <- "No"
chem_CAS <- rbind(chem_CAS_non_toxcast,chem_CAS_toxcast)

files_toxcast <- list.files(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ToxCast"), pattern="*.txt", all.files=FALSE, full.names=FALSE)
files_non_toxcast <- "ECOTOX_nontoxcast.txt"
files <- c(files_non_toxcast,files_toxcast)

chems <- chem_CAS$chnm
chem_CAS$CAS.Number. <- as.numeric(gsub("-","",chem_CAS$CAS))


for (i in 1:length(files)) {
  if(i == 1) {tox_temp <- read.delim(file = file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","non-ToxCast",files[i]),sep="|", stringsAsFactors = FALSE)
  }else {tox_temp <- read.delim(file = file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","ToxCast",files[i]),sep="|", stringsAsFactors = FALSE)
  }
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

remove_references <- c(67566, 168095, 171681, 168095, 168095, 11170, 11628,160420,89736,174456)#, 157699)
# Carbamazepine (157699) moderate outlier: mRNA signals. Removed per Dan and Brett. All others are fairly large outliers.

tox_fw <- tox %>%
  filter(Media.Type == "Fresh water",
         Effect != "Accumulation",
         Exposure.Type %in% exposure.type.keep,
         Conc.1.Type..Standardized..  %in% c("Active ingredient","Total"),
         grepl("mg/L|ug/L",Conc.1.Units..Standardized..),
         !grepl("No significance",Statistical.Significance.),
         !(value < 4.80E-9 & CAS == "1912-24-9"), #Atrazine outlier
         !(value < 0.062 & CAS == "21145-77-7"),  # Tonalide outlier
         !(Reference.Number. %in% remove_references))  #Pyrene study with sediment and pore water

ugL <- grep("ug/L",tox_fw$Conc.1.Units..Standardized..)

if(length(ugL) > 0) {tox_fw$value[ugL] <- tox_fw$value[ugL]/1000} #convert ug/L to mg/L
  
grep("ug/L",tox_fw$Conc.1.Units..Standardized..)
test <- tox[grep("ug/L",tox$Conc.1.Units..Standardized..),]

tox_fw <- tox_fw %>%
  arrange(chnm,value) 

num_chems <- numeric()
chem_index <- numeric()
for (i in 1:length(unique(tox_fw$chnm))) {
  num_chems <- sum(tox_fw$chnm == unique(tox_fw$chnm)[i])
  chem_index <- c(chem_index,1:num_chems)
}
tox_fw$index <- chem_index

All_effects <- unique(tox$Effect)
EffectCategory1 <- c("Reproduction","Mortality","Growth","Development","Population", "Behavior")
EffectCategory2 <- All_effects[!(All_effects %in% EffectCategory1)]
tox_fw$EffectCategory <- ifelse(tox_fw$Effect %in% EffectCategory1,1,2)



benchmark_tab <- tox_fw[,c("CAS.Number.","Chemical.Name","value", "Observed.Duration.Mean..Days..", "Endpoint","Effect","Effect.Measurement","EffectCategory")]
names(benchmark_tab) <- c("CAS.Number.","Chemical.Name","Value", "duration", "Endpoint_type","Effect","Effect.Measurement","groupCol")

benchmark_tab <- benchmark_tab %>%
  mutate(endPoint = ifelse(duration > 4,"Chronic","Acute")) 

benchmark_tab <- left_join(benchmark_tab,chem_CAS[,c("CAS.Number.", "CAS","chnm")]) %>%
  rename(Chemical = chnm)

path_to_data <- Sys.getenv("PASSIVE_PATH")

wb <- loadWorkbook(file.path(path_to_data, "data", "data_for_git_repo","clean", "passive.xlsx"))
#saveWorkbook(wb,file=file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_TEST.xlsx"),overwrite = TRUE)
addWorksheet(wb,sheetName = "Benchmarks")
writeData(wb,sheet = "Benchmarks",x=benchmark_tab)
saveWorkbook(wb,file=file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_all.xlsx"),overwrite = TRUE)
dir.create("R/Analyze/Out", showWarnings = FALSE)
write.csv(tox_fw,"R/Analyze/Out/ECOTOX_combined.csv",row.names = FALSE)
saveRDS(tox_fw,"R/Analyze/Out/ECOTOX_combined.Rds")
