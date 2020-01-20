library(tidyverse)

files <- list.files(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX"), pattern="*.txt", all.files=FALSE, full.names=FALSE)

chems <- sub(x = files,".txt","",)


for (i in 1:length(files)) {
  tox_temp <- read.delim(file = file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX",files[i]),sep="|", stringsAsFactors = FALSE)
  names(tox_temp) <- sub("X.","",names(tox_temp))
  tox_temp$Conc.1.Mean..Standardized.. <- as.numeric(tox_temp$Conc.1.Mean..Standardized..)
  tox_temp$Conc.Min.1..Standardized.. <- as.numeric(tox_temp$Conc.Min.1..Standardized..)
  tox_temp$Observed.Duration.Mean..Days.. <- as.numeric(tox_temp$Observed.Duration.Mean..Days..)
  tox_temp$Chemical.Purity.Mean.... <- as.numeric(tox_temp$Chemical.Purity.Mean....)
  tox_temp$Organism.Age.Mean. <- as.numeric(tox_temp$Organism.Age.Mean.)
  tox_temp$Conc.1.Max..Standardized.. <- as.numeric(tox_temp$Conc.1.Max..Standardized..)
  tox_temp$chems <- chems[i]
  if(i == 1) {tox <- tox_temp 
  }else tox <- bind_rows(tox,tox_temp)
}


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
  arrange(chems,value) 
# %>%
#   group_by(chems) %>%
#   transform(Endpoint_num = 1:length(chems))
# transform(Endpoint_num = 1:length(chems))

num_chems <- numeric()
chem_index <- numeric()
for (i in 1:length(unique(tox_fw$chems))) {
  num_chems <- sum(tox_fw$chems == unique(tox_fw$chems)[i])
  chem_index <- c(chem_index,1:num_chems)
}
tox_fw$index <- chem_index


ggplot(data = tox_fw,aes(x=Effect,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chems)

#Cumulative distribution curve below:
#There do not look to be any anamalously low values, so use the minimum value 
#for each chemical to compare against.
ggplot(data = tox_fw,aes(x=index,y=value)) + 
  geom_point()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chems, scales="free_x")

ggplot(data = tox_fw,aes(x=value,y=index)) + 
  geom_point()+
  scale_x_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chems, scales="free_y")


ggplot(data = tox_fw,aes(x=chems,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90))

# Subset cis and trans isomers and rerun graphs:
#Cumulative distribution curve below:
#There do not look to be any anamalously low values, so use the minimum value 
#for each chemical to compare against.
tox_isomers <- tox_fw[grep("cis|trans",tox_fw$chems,ignore.case = TRUE),]
unique(tox_isomers$chems)

tox_isomers$chems <- factor(tox_isomers$chems,
                            levels = c("cis-chlordane","trans-chlordane",
                                       "trans-nonachlor",
                                       "cis-permethrin","trans-permethrin"))
ggplot(data = tox_isomers,aes(x=index,y=value)) + 
  geom_point()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chems, scales="free_x")

ggplot(data = tox_isomers,aes(x=value,y=index)) + 
  geom_point()+
  scale_x_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chems, scales="free_y")


ggplot(data = tox_isomers,aes(x=chems,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90))


c("Fish","Algae","Amphibians","Crustaceans","Flowers","insects/Spiders","Invertebrates","Molluscs")
