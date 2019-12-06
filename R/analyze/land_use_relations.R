#Explore land cover relations with different contaminant classes

library(toxEval)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)

# -Classes that are likely to have EAR > 0.001
#   -Endpoints that are important for these classes
# -Full mixtures
#   -Endpoints that are important for full mixtures
#   -Chemicals/classes that play prominent roles in the EARs from mixtures

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




########################################################################################
# 1). Determine EARs summation for this analysis using by endpoints within chemical classes 
# 2). Find max of those EAR summations

# Determine EAR_sum by site/class/endpoint
class_EP_summary <- chemical_summary %>%
  group_by(site,date,Class, endPoint) %>%
  summarize(EAR_sum = sum(EAR)) %>%
  group_by(site,Class) 


# Determine max EAR_sum (summed by endpoint/class/site)
class_EP_max <- class_EP_summary %>% 
  summarize(EAR_max = max(EAR_sum),
            EP_max = endPoint[which.max(EAR_sum)])

# Add land use
lu_columns <- c("STAID", "Urban.......6","Parking.Lot....","Agriculture.......7","Crops....","Water.......14","Wetlands....","Population.Density..people.km2.")
names(lu_columns) <- c("site","Urban","Parking_lot","Agriculture","Crops","Water","Wetlands","Population_density")
class_EP_max <- left_join(class_EP_max,df_lu[,lu_columns],by=c("site"="STAID"))
names(class_EP_max)[c(1,5:11)] <- names(lu_columns)


lu_eval <- "Wetlands"
class_EP_max$gtthresh <- ifelse(class_EP_max$EAR_max >= 0.001,"EAR >= 0.001","EAR < 0.001")
gg.summary <- group_by(class_EP_max, Class,gtthresh) %>% summarise(length=length(EP_max))  
p_urban <- ggplot(data=class_EP_max,aes_string(x="gtthresh",y=lu_eval)) +
  geom_boxplot() + 
  #  scale_y_continuous(trans='log10') +
  facet_wrap(~Class,nrow = 3) +
  geom_text(data = gg.summary,
            aes(gtthresh, Inf, label = length), size=3,vjust = 3,) + 
  ggtitle(lu_eval) + 
  ylab(paste("%",lu_eval,"Land Cover")) + 
  #  ylim(c(0,120)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,150)) +
  stat_compare_means()

p_urban

lu_eval <- "Crops"
class_EP_max$gtthresh <- ifelse(class_EP_max$EAR_max >= 0.001,"EAR >= 0.001","EAR < 0.001")
gg.summary <- group_by(class_EP_max, Class,gtthresh) %>% summarise(length=length(EP_max))  
p_crops <- ggplot(data=class_EP_max,aes_string(x="gtthresh",y=lu_eval)) +
  geom_boxplot() +
  #  scale_y_continuous(trans='log10') +
  facet_wrap(~Class,nrow = 3) +
  geom_text(data = gg.summary,
            aes(gtthresh, Inf, label = length), size=3,vjust = 3,) +
  ggtitle(lu_eval) +
  ylab(paste("%",lu_eval,"Land Cover")) +
  #  ylim(c(0,120)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,150)) +
  stat_compare_means()

p_crops

lu_eval <- "Agriculture"
class_EP_max$gtthresh <- ifelse(class_EP_max$EAR_max >= 0.001,"EAR >= 0.001","EAR < 0.001")
gg.summary <- group_by(class_EP_max, Class,gtthresh) %>% summarise(length=length(EP_max))  
p_Ag <- ggplot(data=class_EP_max,aes_string(x="gtthresh",y=lu_eval)) +
  geom_boxplot() + 
  #  scale_y_continuous(trans='log10') +
  facet_wrap(~Class,nrow = 3) +
  geom_text(data = gg.summary,
            aes(gtthresh, Inf, label = length), size=3,vjust = 3,) + 
  ggtitle(lu_eval) + 
  ylab(paste("%",lu_eval,"Land Cover")) + 
  #  ylim(c(0,120)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,150)) +
  stat_compare_means()

p_Ag


class_EP_max$Developed <- class_EP_max$Urban + class_EP_max$Agriculture
lu_eval <- "Developed"
class_EP_max$gtthresh <- ifelse(class_EP_max$EAR_max >= 0.001,"EAR >= 0.001","EAR < 0.001")
gg.summary <- group_by(class_EP_max, Class,gtthresh) %>% summarise(length=length(EP_max))  
p_Dev <- ggplot(data=class_EP_max,aes_string(x="gtthresh",y=lu_eval)) +
  geom_boxplot() + 
  #  scale_y_continuous(trans='log10') +
  facet_wrap(~Class,nrow = 3) +
  geom_text(data = gg.summary,
            aes(gtthresh, Inf, label = length), size=3,vjust = 3,) + 
  ggtitle(lu_eval) + 
  ylab(paste("%",lu_eval,"Land Cover")) + 
  #  ylim(c(0,120)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,150)) +
  stat_compare_means()

p_Dev



####################################################################################
## Explore the endpoints 


thresh <- 0.001
class_EP_summary_thresh <- filter(class_EP_summary,EAR_sum > 0.001)
length(unique(test$endPoint))

class_list <- unique(class_EP_summary_thresh$Class)

filenm <- "endpoints_by_class.pdf"
pdf()
for (i in 1:length(class_list)) {
  
  plot_df <- filter(class_EP_summary_thresh,Class == class_list[i])
    
  p <- ggplot(data = plot_df,aes(x=EAR_sum,y=endPoint)) +
    geom_boxplot() +
    facet_wrap(~Class) +
    scale_x_continuous(trans='log10')
  print(p)
}
dev.off()
shell.exec(filenm)


####################################
## LU vs EAR by Class Regressions ##


df <- class_EP_max %>% 
  filter(Class == "PAHs")


summary(lm(EAR_max~Urban + Agriculture + Water + Wetlands, data=df))
summary(lm(EAR_max~Urban + Agriculture, data=df))
summary(lm(EAR_max~Urban , data=df))
summary(lm(EAR_max~Urban + Crops, data=df))
summary(lm(EAR_max~ Crops, data=df))
summary(lm(EAR_max~Crops + Water + Wetlands, data=df))
summary(lm(EAR_max~Urban + Agriculture, data=df))
summary(lm(EAR_max~Population_density, data=df))
summary(lm(EAR_max~Parking_lot, data=df))

lu_variables <- names(lu_columns[-1])
form <- as.formula(paste("EAR_max ~ ",paste(lu_variables,collapse = " + ")))

classes_reg <- c("Herbicide","Fire retardant", "Pharmaceuticals","PAHs")
for (i in 1:length(classes_reg)){
  df <- class_EP_max %>% 
    filter(Class == classes_reg[i])
  m <- lm(EAR_max~1,data=df)
  
  m_step <- step(m,scope = form,data = df,k = log(dim(df)[1]))
  




  # Graph EAR summaries vs land use
# Compute EAR summaries by Class regardless of endpoint

# class_summary <- chemical_summary %>% group_by(site,Class) %>%
#   summarize(EAR_sum = sum(EAR))
# 
# lu_columns <- c("STAID", "Urban.......6","Parking.Lot....","Agriculture.......7")
# names(lu_columns) <- c("site","Urban","Parking_lot","Agriculture")
# 
# class_summary <- left_join(class_summary,df_lu[,lu_columns],by=c("site"="STAID"))
# names(class_summary)[c(1,4,5,6)] <- names(lu_columns)
# 
# unique(class_summary$Class)
# 
# ggplot(data=class_summary,aes(x=Parking_lot,y=EAR_sum)) +
#   geom_point() + 
#   scale_y_continuous(trans='log10') +
#   facet_wrap(~Class)
# 
# plot(class_summary$Urban,class_summary$Parking.Lot.....y)
# 
# 
# class_summary$urban_cat_15 <- ifelse(class_summary$Urban > 15,"high_urban","low_urban")
# ggplot(data=class_summary,aes(x=urban_cat_15,y=EAR_sum)) +
#   geom_boxplot() + 
#   scale_y_continuous(trans='log10') +
#   facet_wrap(~Class)
# 
# gg.summary <- group_by(class_summary, Class,gtthresh) %>% summarise(length=length(EAR_sum))  
# class_summary$gtthresh <- ifelse(class_summary$EAR_sum > 0.001,"exceedance","non-exceedance")
# ggplot(data=class_summary,aes(x=gtthresh,y=Urban)) +
#   geom_boxplot() + 
# #  scale_y_continuous(trans='log10') +
#   facet_wrap(~Class)+
#   geom_text(data = gg.summary,
#             aes(gtthresh, Inf, label = length), size=3,vjust = 0.5)
# 
# class_summary$gtthresh <- ifelse(class_summary$EAR_sum > 0.001,"exceedance","non-exceedance")
# p <- ggplot(data=class_summary,aes(x=gtthresh,y=Agriculture)) +
#   geom_boxplot() + 
#   #  scale_y_continuous(trans='log10') +
#   facet_wrap(~Class)
# 
# 
# gg.summary <- group_by(class_summary, Class,gtthresh) %>% summarise(length=length(EAR_sum))  
# class_summary$gtthresh <- ifelse(class_summary$EAR_sum > 0.001,"exceedance","non-exceedance")
# p <- ggplot(data=class_summary,aes(x=gtthresh,y=Agriculture)) +
#   geom_boxplot() + 
#   #  scale_y_continuous(trans='log10') +
#   facet_wrap(~Class) +
#   geom_text(data = gg.summary,
#             aes(gtthresh, Inf, label = length), size=3,vjust = 1.1,)
# p
