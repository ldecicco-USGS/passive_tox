#Vetting ECOTOX endpoints
library(tidyverse)

#Subset lowest values for vetting of priorities (> thresh at x% of sites)
site_detections <- read.csv("R/Analyze/Out/ECOTOX_site_threshold_exceedances_questionable.csv",stringsAsFactors = FALSE)

tox_validation <- tox_fw[0:0,]
for(i in 1:dim(site_detections)[1] ) {
  CASnum <- site_detections[i,"CAS"]
  conc <- site_detections[i,"Check.concentration"]
  tox_validation <- suppressMessages( bind_rows(tox_validation,filter(tox_fw,(CAS == CASnum & value < conc))))
}

tox_validation <- tox_validation %>%
  arrange(chnm,value)

#tox_validation <- tox_validation[,c("CAS.Number.","Chemical.Name","Species.Scientific.Name.","Species.Common.Name","Species.Group"

unique(tox_validation$CAS)
unique(tox_validation$Effect)
table(tox_validation$Statistical.Significance.)

table(tox_validation$chnm)
write.csv(tox_validation,"R/Analyze/Out/ECOTOX_results_questionable.csv")
