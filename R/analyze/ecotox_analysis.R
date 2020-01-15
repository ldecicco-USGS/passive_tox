library(tidyverse)

files <- list.files(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX"), pattern="*.txt", all.files=FALSE, full.names=FALSE)

chems <- sub(x = files,".txt","",)

for (i in 1:length(files)) {
  tox_temp <- read.delim(file = file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX",files[i]),sep="|", stringsAsFactors = FALSE)
  names(tox_temp) <- sub("X.","",names(tox_temp))
  tox_temp$Conc.1.Mean..Standardized.. <- as.numeric(tox_temp$Conc.1.Mean..Standardized..)
  tox_temp$Conc.1.Min..Standardized.. <- as.numeric(tox_temp$Conc.1.Min..Standardized..)
  tox_temp$chems <- chems[i]
  if(i == 1) {tox <- tox_temp 
  }else tox <- bind_rows(tox,tox_temp)
}

tox$Conc.1.Mean..Standardized.. <- as.numeric(tox$Conc.1.Mean..Standardized..)
tox$Conc.Min.1..Standardized.. <- as.numeric(tox$Conc.Min.1..Standardized..)
tox <- transform(tox,value = pmin(Conc.1.Mean..Standardized..,Conc.Min.1..Standardized..,na.rm=TRUE))
sum(is.na(tox$value))

test <- ifelse(is.na(tox$Conc.Min.1..Standardized..),tox$Conc.1.Mean..Standardized..,tox$Conc.Min.1..Standardized..)
sum(is.na(test))

tox_fw <- tox %>%
  filter(Media.Type == "Fresh water",
         Conc.1.Type..Standardized.. == "Active ingredient")

boxplot(value~Media.Type,log="y",las=2,data=tox)

ggplot(data = tox_fw,aes(x=Effect,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90))

ggplot(data = tox_fw,aes(x=Chemical.Name,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90))


c("Fish","Algae","Amphibians","Crustaceans","Flowers","insects/Spiders","Invertebrates","Molluscs")
