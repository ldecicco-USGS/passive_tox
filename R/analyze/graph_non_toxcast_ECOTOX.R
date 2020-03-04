# Graph non-toxcast ECOTOX analysis results

tox_fw <- readRDS("R/Analyze/Out/ECOTOX_filtered_non_toxcast.Rds")


chem_CAS <- read.csv(file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","non-ToxCast","chem_info_non_toxcast.csv"),stringsAsFactors = FALSE)

num_cols <- 3
num_rows <- 4
pages <-  ceiling(length(unique(tox_fw$CAS.Number.))/(num_cols*num_rows))


pdf("R/Analyze/Plots/non_ToxCast_chems_boxplots.pdf")
for(i in 1:pages) {
  p <- ggplot(data = tox_fw,aes(x=Effect,y=value)) + 
    geom_boxplot()+
    scale_y_continuous(trans='log10') + 
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap_paginate(~chnm, ncol = num_cols, nrow = num_rows,page = i)
  
  print(p)
}
dev.off()


#Cumulative distribution curve below:
#There do not look to be any anamalously low values, so use the minimum value 
#for each chemical to compare against.
pdf("R/Analyze/Plots/non_ToxCast_chems_CDF.pdf")
for(i in 1:pages) {
  CDF <- ggplot(data = tox_fw,aes(y=index,x=value)) + 
    geom_point()+
    scale_x_continuous(trans='log10') + 
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap_paginate(~chnm, ncol=num_cols,nrow=num_rows,page = i,scales="free")
  
  
  print(CDF)
}
dev.off()



#Cumulative distribution curve below:
#There do not look to be any anamalously low values, so use the minimum value 
#for each chemical to compare against.
CDF <- ggplot(data = tox_fw,aes(x=index,y=value)) + 
  geom_point()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~chnm, scales="free_x")
CDF


ggplot(data = tox_fw,aes(x=chnm,y=value)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10') + 
  theme(axis.text.x = element_text(angle = 90))

