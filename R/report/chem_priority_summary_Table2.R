

library(toxEval)
library(tidyverse)

get_table_2 <- function(){
  path_to_data <- Sys.getenv("PASSIVE_PATH")
  
  source("read_chemicalSummary.R")
  chem_data <- tox_list$chem_data
  chem_info <- tox_list$chem_info
  
  #Read individual priority chems and mixture priority chems
  priority_chems <- readRDS("R/analyze/out/priority_chem_EAR_TQ.rds")
  mixtures <- readRDS("R/mixtures/out/mixtures_table.rds")
  names(mixtures)[8]

  mix_chems <- data.frame()
  for(i in 1:dim(mixtures)[1]) {
    row_chems <- unlist(mixtures$Chemicals[i])
    row_CAS <- unlist(mixtures$CASs[i])
    row_sites <- unlist(mixtures[i,8])
    mix_chems <- rbind(mix_chems,
                       data.frame(mix_chemicals = row_chems,
                                  CAS = row_CAS,
                                  Mixtures = rep(row_sites,length(row_CAS)),
                                  stringsAsFactors = FALSE))
  }
  
  mix_chems <- mix_chems %>%
    group_by(mix_chemicals,CAS) %>%
    summarize(Mixtures = max(Mixtures))
  
  chem_priorities <- full_join(mix_chems,priority_chems)
  
  chem_priorities$Chemicals <- ifelse(is.na(chem_priorities$Mixtures),chem_priorities$chnm,as.character(chem_priorities$mix_chemicals))
  
  names(chem_priorities)
  
  table_2 <- chem_priorities %>%
    filter(mix_chemicals > 0 | max_exceed >= 0.1) %>%
    ungroup()%>%
    select(Class,Chemicals,sites_monitored,ToxCast,Mixtures,ECOTOX_group_1,ECOTOX_group_2) %>%
    arrange(Class,Chemicals)
  
  #Populate mixtures columns for individual chems that have missing info
  add_nums_to_mixtures <- which(is.na(table_2$Mixtures) & !is.na(table_2$ToxCast))
  table_2[add_nums_to_mixtures,"Mixtures"] <- table_2[add_nums_to_mixtures,"ToxCast"]
  
  table_2[, 4:7] <- sapply(table_2[, 4:7], as.character)
  for (i in 4:7) {
    table_2[which(is.na(table_2[,i])),i] <- "--"
  }
  
  #add empty column for spacer
  table_2$empty <- ""
  table_2 <- table_2[,c(1:5,8,6:7)]
  
  return(table_2)
}



