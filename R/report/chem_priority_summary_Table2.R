library(toxEval)
library(tidyverse)

path_to_data <- Sys.getenv("PASSIVE_PATH")

source("read_chemicalSummary.R")
chem_data <- tox_list$chem_data
chem_info <- tox_list$chem_info

priority_chems <- readRDS("R/analyze/out/priority_chem_EAR_TQ.rds")
mixture_chems <- read.csv(file.path(path_to_data,"/Tables/mixture_chems.csv"),stringsAsFactors = FALSE)
mixture_chems$mixture <- "X"
mixture_chems <- mixture_chems[,-1]
names(mixture_chems)[1] <- "Chemicals"
names(priority_chems)

#keep_cols <- c("Class","chnm","sites_monitored","ToxCast","ECOTOX_group_1","ECOTOX_group_2")

#table2 <- priority_chems[,keep_cols]

mixture_chems[which(!(mixture_chems$Chemicals %in% chem_info$chnm)),"Chemicals"]

mixture_chems[grep("citrate",mixture_chems$Chemicals,ignore.case = TRUE),"Chemicals"]
chem_info[grep("citrate",chem_info$chnm,ignore.case = TRUE),"chnm"]







priority_chems$Chemicals <- priority_chems$chnm


chem.chnm <- c("Diethylhexylphthalate (DEHP)", "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)", 
       "Tris(2-butoxyethyl)phosphate (TBEP)", "Tris(1-chloro-2-propyl)phosphate (TCPP)",
       "N,N-diethyltoluamide (DEET)")
chem.Chemicals <- c("DEHP", "TDCPP", 
                    "TBEP", "TCPP", "DEET")
       
priority_chems$Chemicals <- priority_chems$chnm
priority_chems$Chemicals[which(priority_chems$chnm %in% chem.chnm)] <- chem.Chemicals


mix.chnm <- c("TDCPP", "TBEP", "TCPP")
mix.Chemicals <- c("Tris(1,3-dichloro-2-propyl)phosphate", 
                    "Tris(2-butoxyethyl)phosphate", "Tris(1-chloro-2-propyl)phosphate")
mixture_chems[which(mixture_chems$Chemicals %in% mix.chnm),"Chemicals"] <- mix.Chemicals


priority_chems_joined <- full_join(mixture_chems,priority_chems)

priority_chems_joined$ToxCast
