#NOTE: Add path to path_to_file!!!
library(toxEval)
library(dplyr)
library(readxl)
library(ggplot2)

df_lu <- read_xlsx(path = file.path("data","raw","GLRItox_summary.xlsx"),sheet = 1,skip=1)
names(df_lu) <- make.names(names(df_lu))

path_to_file <- 'passive.xlsx' 
tox_list <- create_toxEval(file.path("data","clean",path_to_file))
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC,
                    flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

class_summary <- chemical_summary %>% group_by(site,Class) %>%
  summarize(EAR_sum = sum(EAR))

lu_columns <- c("STAID", "Urban.......6","Parking.Lot....","Agriculture.......7")
names(lu_columns) <- c("site","Urban","Parking_lot","Agriculture")

class_summary <- left_join(class_summary,df_lu[,lu_columns],by=c("site"="STAID"))
names(class_summary)[c(1,4,5,6)] <- names(lu_columns)

unique(class_summary$Class)

ggplot(data=class_summary,aes(x=Parking_lot,y=EAR_sum)) +
  geom_point() + 
  scale_y_continuous(trans='log10') +
  facet_wrap(~Class)

plot(class_summary$Urban,class_summary$Parking.Lot.....y)


class_summary$urban_cat_15 <- ifelse(class_summary$Urban > 15,"high_urban","low_urban")
ggplot(data=class_summary,aes(x=urban_cat_15,y=EAR_sum)) +
  geom_boxplot() + 
  scale_y_continuous(trans='log10') +
  facet_wrap(~Class)

gg.summary <- group_by(class_summary, Class,gtthresh) %>% summarise(length=length(EAR_sum))  
class_summary$gtthresh <- ifelse(class_summary$EAR_sum > 0.001,"exceedance","non-exceedance")
ggplot(data=class_summary,aes(x=gtthresh,y=Urban)) +
  geom_boxplot() + 
#  scale_y_continuous(trans='log10') +
  facet_wrap(~Class)+
  geom_text(data = gg.summary,
            aes(gtthresh, Inf, label = length), size=3,vjust = 0.5)

class_summary$gtthresh <- ifelse(class_summary$EAR_sum > 0.001,"exceedance","non-exceedance")
p <- ggplot(data=class_summary,aes(x=gtthresh,y=Agriculture)) +
  geom_boxplot() + 
  #  scale_y_continuous(trans='log10') +
  facet_wrap(~Class)


gg.summary <- group_by(class_summary, Class,gtthresh) %>% summarise(length=length(EAR_sum))  
class_summary$gtthresh <- ifelse(class_summary$EAR_sum > 0.001,"exceedance","non-exceedance")
p <- ggplot(data=class_summary,aes(x=gtthresh,y=Agriculture)) +
  geom_boxplot() + 
  #  scale_y_continuous(trans='log10') +
  facet_wrap(~Class) +
  geom_text(data = gg.summary,
            aes(gtthresh, Inf, label = length), size=3,vjust = 1.1,)
p
