library(tidyverse)
library(toxEval)
library(ggforce)

source(file = "read_chemicalSummary.R")
source(file = "R/mixtures/mix_script.R")

#Get some info on ToxCast assays
# Original assays:
ToxCast_ACC <- ToxCast_ACC

# tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
#                                      "data","data_for_git_repo","clean",
#                                      "passive.xlsx"))

casn <- tox_list$chem_info$CAS
chemicalSummary <- chemicalSummary %>%
  filter(EAR > 0)
endpoints <- unique(chemicalSummary$endPoint)

ACC <- filter(ToxCast_ACC, CAS %in% casn,endPoint %in% endpoints)
  


ACC_min <- ACC %>%
  group_by(CAS) %>%
  summarize(ACC_min = min(ACC)) %>%
  left_join(tox_list$chem_info[,c("CAS","chnm","Class")]) %>%
  arrange(chnm)

ACC <- left_join(ACC,tox_list$chem_info)
ACC$chnm <- factor(ACC$chnm,levels = ACC_min$chnm)

EP_group <- character()
for(i in 1:dim(ACC)[1]){
  EP_group <- c(EP_group,strsplit(ACC$endPoint[i],split = "_")[[1]][1])
}
ACC$EP_group <- EP_group
ACC <- arrange(ACC,chnm,ACC)



# Plot 12 chems per page as scatter of ACC values
nchems <- length(unique(ACC$chnm))
pages <- ceiling(nchems/12)

pdf(file="./R/analyze/Plots/ACC_scatter_ecdf.pdf")
for(i in 1:pages) {
  p<- ggplot(ACC,aes(x=ACC)) + 
    stat_ecdf(geom = "point", size = 1.5) +
    facet_wrap_paginate(.~chnm,scales = "free",nrow=4,ncol=3,page = i) 
  print(p)
}    

dev.off()
