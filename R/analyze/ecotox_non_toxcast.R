library(tidyverse)
library(toxEval)
library(readxl)


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

tox <- transform(tox,value = pmin(conc_mean,conc_min,conc_max,na.rm=TRUE))
sum(is.na(tox$value))

tox_fw <- tox %>%
  filter(Media.Type == "Fresh water",
         Conc.1.Type..Standardized.. == "Active ingredient",
         Effect != "Accumulation",
         Exposure.Type != "Intraperitoneal",
         Exposure.Type != "Food",
         Exposure.Type != "Injection, unspecified",
         Exposure.Type != "Intramuscular",
         Media.Type != "Salt water",
         Conc.1.Units..Standardized.. %in% c("AI mg/L", "ml/L"))

# Statistical.Significance != "No significance",
         # Statistical.Significance != "Not significant at all concentrations")

boxplot(value~Media.Type,log="y",las=2,data=tox)

test <- tox[which(!(tox$value > 0)),]

tox_fw <- tox_fw %>%
  arrange(chnm,value) 
# %>%
#   group_by(chnm) %>%
#   transform(Endpoint_num = 1:length(chnm))
# transform(Endpoint_num = 1:length(chnm))

num_chems <- numeric()
chem_index <- numeric()
for (i in 1:length(unique(tox_fw$chnm))) {
  num_chems <- sum(tox_fw$chnm == unique(tox_fw$chnm)[i])
  chem_index <- c(chem_index,1:num_chems)
}
tox_fw$index <- chem_index

# num_chems <- numeric()
# chem_index <- numeric()
# chem_CAS$CAS.Number. <- gsub("-","",x = chem_CAS$CAS)
# 
# 
# for (i in 1:length(unique(tox_fw$Chemical.Name))) {
#   num_chems <- sum(tox_fw$Chemical.Name == unique(tox_fw$Chemical.Name)[i])
#   chem_index <- c(chem_index,1:num_chems)
# }
# tox_fw$index <- chem_index

ggplot(data = tox_fw,aes(x=Effect,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chnm)

#Cumulative distribution curve below:
#There do not look to be any anamalously low values, so use the minimum value 
#for each chemical to compare against.
CDF <- ggplot(data = tox_fw,aes(x=index,y=value)) + 
  geom_point()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chnm, scales="free_x")
CDF

# CDF2 <- CDF +
#   geom_point(data=ACC)+
#   theme(axis.text.x = element_text(angle = 90)) +
#   facet_wrap(~chnm, scales="free_x")
# 
# CDF2
# 
ggplot(data = tox_fw,aes(x=value,y=index)) + 
  geom_point()+
  scale_x_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chnm, scales="free_y")


ggplot(data = tox_fw,aes(x=chnm,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90))


#Determine stats for each chem
tox_stats <- tox_fw %>%
  group_by(chnm,CAS.Number.,Chemical.Name) %>%
  summarize(min_endpoint = min(value),
            median_endpoint = median(value),
            num_endpoints = length(unique(value)))

c("Fish","Algae","Amphibians","Crustaceans","Flowers","insects/Spiders","Invertebrates","Molluscs")
