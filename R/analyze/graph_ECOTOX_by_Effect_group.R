library(tidyverse)
library(toxEval)
library(ggforce)

#source(file = "read_chemicalSummary.R")
#source(file = "R/mixtures/mix_script.R")


tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","toxEval input file",
                                     "passive_benchmarks_all.xlsx"))
benchmarks <- tox_list$benchmarks

endpoints <- unique(benchmarks$endPoint)
casn <- unique(benchmarks$CAS.Number.)

ACC <- filter(benchmarks,CAS.Number. %in% casn,endPoint %in% endpoints)
ACC$value <- ACC$ACC_value
ACC$ACC <- ACC$ACC_value

ACC_min <- ACC %>%
  group_by(CAS) %>%
  summarize(ACC_min = min(value)) %>%
  left_join(tox_list$chem_info[,c("CAS","chnm","Class")]) %>%
  arrange(chnm)


ACC$EP_group <- ACC$Effect
tier1 <- c("Reproduction","Mortality","Growth","Development","Population", "Behavior")

ACC$tier <- ifelse(ACC$EP_group %in% tier1,"Tier 1","Tier 2")

# Plot 12 chems per page as scatter of ACC values
nchems <- length(unique(ACC$chnm))
pages <- ceiling(nchems/12)
pdf(file="./R/analyze/Plots/ECOTOX_scatter.pdf")
for(i in 1:pages) {
  p<- ggplot(ACC,aes(x=ACC,y=EP_group,colour = tier)) + 
    geom_point() + 
    facet_wrap_paginate(.~chnm,scales = "free",nrow=4,ncol=3,page = i) +
    theme(axis.text.y=element_text(angle=0, size=5, vjust=0.5),
          axis.text.x=element_text(angle=90, size=5, vjust=0.5)) + 
    scale_x_log10()
  print(p)
}    

dev.off()
