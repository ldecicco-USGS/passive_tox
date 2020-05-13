# Combine ECOTOX analysis for toxcast and non toxcast chems

library(openxlsx)
library(readxl)

ecotox_toxcast <- read_rds("R/Analyze/Out/ECOTOX_filtered_toxcast.Rds")
ecotox_non_toxcast <- read_rds("R/Analyze/Out/ECOTOX_filtered_non_toxcast.Rds")
ecotox_toxcast <- read_rds("R/Analyze/Out/ECOTOX_filtered_toxcast.Rds")
ecotox_non_toxcast <- read_rds("R/Analyze/Out/ECOTOX_filtered_non_toxcast.Rds")

#
# Combine the toxEval benchmarks tabs and write comprehensive toxEval input file.
wb_toxcast <- loadWorkbook(file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx"))
wb_non_toxcast <- loadWorkbook(file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx"))
bm_toxcast <- read_xlsx(file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx"),sheet = "Benchmarks")
bm_toxcast$ToxCast_inclusion <- "Yes"
bm_non_toxcast <- read_xlsx(file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx"),sheet = "Benchmarks")
bm_non_toxcast$ToxCast_inclusion <- "No"
bm <- full_join(bm_toxcast,bm_non_toxcast)
wb <- wb_toxcast
writeData(wb,sheet = "Benchmarks",x=bm)
saveWorkbook(wb,file=file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks.xlsx"),overwrite = TRUE)


#CombineECOTOX stat summaries for SI table
tox_stats_toxcast <- readRDS(file = "R/Analyze/Out/Tox_endpoint_stats_toxcast.rds")
tox_stats_non_toxcast <- readRDS(file = "R/Analyze/Out/Tox_endpoint_stats_non_toxcast.rds")
tox_stats_toxcast$ToxCast_inclusion <- "Yes"
tox_stats_non_toxcast$ToxCast_inclusion <- "No"
tox_stats <- full_join(tox_stats_toxcast,tox_stats_non_toxcast)
write.csv(tox_stats,file = file.path(Sys.getenv("PASSIVE_PATH"),"Supplemental","SI_table_ECOTOX_stats.csv"),row.names = FALSE)


# Combine full ECOTOX downloads and write to file for SI table
logical_cols <- grep("logical",sapply(ecotox_toxcast,class))
logical_cols <- c(logical_cols,grep("logical",sapply(ecotox_non_toxcast,class)))
logical_cols <- unique(logical_cols)

for (i in logical_cols) {
  ecotox_toxcast[,i] <- as.character(ecotox_toxcast[,i])
  ecotox_non_toxcast[,i] <- as.character(ecotox_non_toxcast[,i])
}

ecotox <- full_join(ecotox_toxcast,ecotox_non_toxcast)
suppressMessages(
include <- ecotox %>%
  select(Class,CAS,chnm,Species.Scientific.Name., Species.Common.Name, Species.Group, 
         Organism.Lifestage, Organism.Age.Mean.Op., Organism.Age.Mean., Organism.Age.Min.Op.,
         Organism.Age.Min., Organism.Age.Max.Op., Organism.Age.Max., Age.Units, Exposure.Type, 
         Media.Type, Test.Location, Number.of.Doses, 
         Conc.1.Type..Standardized.., Conc.1.Mean.Op..Standardized.., Conc.1.Mean..Standardized.., 
         Conc.1.Min.Op..Standardized.., Conc.Min.1..Standardized.., Conc.1.Max.Op..Standardized.., 
         Conc.1.Max..Standardized.., Conc.1.Units..Standardized.., Conc.2.Type..Standardized.., 
         Conc.2.Mean.Op..Standardized.., Conc.2.Mean..Standardized.., Conc.2.Min.Op..Standardized.., 
         Conc.Min.2..Standardized.., Conc.2.Max.Op..Standardized.., Conc.2.Max..Standardized.., 
         Conc.2.Units..Standardized.., Conc.3.Type..Standardized.., Conc.3.Mean.Op..Standardized.., 
         Conc.3.Mean..Standardized.., Conc.3.Min.Op..Standardized.., Conc.Min.3..Standardized.., 
         Conc.3.Max.Op..Standardized.., Conc.3.Max..Standardized.., Conc.3.Units..Standardized., 
         Effect, Effect.Measurement, Endpoint, Response.Site, 
         Statistical.Significance., Significance.Level.Mean.Op., Significance.Level.Mean., 
         Significance.Level.Min.Op., Significance.Level.Min., Significance.Level.Max.Op., Significance.Level.Max, 
         Observed.Duration.Mean.Op..Days.., Observed.Duration.Mean..Days.., Observed.Duration.Min.Op..Days.., 
         Observed.Duration.Min..Days.., Observed.Duration.Max.Op..Days.., Observed.Duration.Max..Days.., 
         Observed.Duration.Units..Days., Author., Reference.Number., Title., Source., 
         Publication.Year, value)
)

write.csv(ecotox,file = file.path("R","Analyze","Out","ECOTOX_combined.csv"),row.names = FALSE)
saveRDS(ecotox,file = file.path("R","Analyze","Out","ECOTOX_combined.rds"))
write.csv(include,file.path(Sys.getenv("PASSIVE_PATH"),"Supplemental","SI_table_ECOTOX.csv"),row.names = FALSE)
