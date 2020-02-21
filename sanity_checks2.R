library(tidyverse)
library(toxEval)
library(openxlsx)

source(file = "read_chemicalSummary.R")
tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))

path_to_data <- Sys.getenv("PASSIVE_PATH")

wb <- loadWorkbook(file.path(path_to_data, "Supplemental", "Table S2 - Chemical list+conc ranges.xlsx"))

df <- readWorkbook(wb, startRow = 5, check.names = TRUE)

df_joined <- df %>% 
  select(Class, Analyte, CAS) %>% 
  filter(Analyte != "") %>% 
  full_join(tox_list$chem_info, by = "CAS") %>% 
  select(CAS, Class.x, Class.y)

changes = data.frame(
  chems = c("Dibenzo[a,h]anthracene", "Bupropion", "Citalopram",
            "Duloxetine", "Methadone", "	Propranolol",
            "Sertraline", "tramadol", "	Tris(1-chloro-2-propyl)phosphate (TCPP)"),
  Dave = c("53-07-3", "34911-55-2", "59729-33-8",
                    "116539-59-4", "76-99-3", "525-66-6",
                    "79617-96-2", "27203-92-5", "26248-87-3"),
  Laura = c("53-70-3", "34841-39-9", "219861-08-2",
                 "136434-34-9", "1095-90-5", "318-98-9",
                 "79559-97-0", "36282-47-0", "13674-84-5"),
  stringsAsFactors = FALSE)

# write.csv(changes, "change_cas.csv", row.names = FALSE)

df <- df %>% 
  filter(!is.na(Analyte))

new_CAS <- df %>% 
  select(CAS) %>% 
  left_join(changes, by=c("CAS"="Dave")) %>% 
  mutate(CAS = ifelse(!is.na(Laura), Laura, CAS)) %>% 
  filter(!is.na(CAS)) 

orig_cas <- df$CAS

df$CAS <- new_CAS$CAS

writeData(wb, sheet = "Table S2",
          startRow = 6, startCol = 3, 
          x = df$CAS)

changedCAS <- createStyle(fgFill = "steelblue2")

addStyle(wb, sheet = "Table S2",
         style = changedCAS,
         cols = 3, gridExpand = FALSE, 
         rows = 5 + which(!orig_cas == new_CAS))


ALL_TOX_DATA <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                  "data","data_for_git_repo","raw",
                                  "all_tox.rds"))

chem_info <- tox_list$chem_info

ALL_TOX_DATA_in_study <-  ALL_TOX_DATA %>% 
  select(CAS = casn, endPoint=aenm, modl_acc, flags, hitc) %>% 
  filter(CAS %in% chem_info$CAS) %>% 
  group_by(CAS) %>% 
  summarize(`Total ToxCast assays` = length(unique(endPoint)),
            `Assays with hits` = length(unique(endPoint[hitc == 1])))

assays_left <- chemicalSummary %>% 
  select(CAS, endPoint) %>% 
  distinct() %>% 
  group_by(CAS) %>% 
  summarize(`Assays in study` = length(unique(endPoint)))

chem_info <- chem_info %>% 
  left_join(ALL_TOX_DATA_in_study, by="CAS") %>% 
  left_join(assays_left, by="CAS")

df_tox <- df %>% 
  select(CAS, Analyte) %>% 
  left_join(chem_info, by="CAS") %>% 
  select(-CAS, -Analyte, -chnm, -Class)

writeData(wb, sheet = "Table S2",
          startRow = 5, startCol = 24, 
          x = df_tox)

chem_stats <- tox_list$chem_data %>% 
  mutate(Value = 1000 * Value) %>% 
  group_by(CAS, `Sample Date`) %>% 
  summarise(min = min(Value[Value > 0]),
            max = max(Value),
            mean = mean(Value[Value > 0]),
            median = median(Value[Value > 0]),
            n_dets = sum(is.na(comment)),
            samples = n()) %>% 
  ungroup() %>% 
  mutate(max = ifelse(max == 0, "ND", 
                      format(round(max, digits = 1), scientific = FALSE)),
         min = ifelse(min == Inf, "ND",
                      format(round(min, digits = 1), scientific = FALSE)),
         mean = ifelse(is.nan(mean), "ND",
                      format(round(mean, digits = 1), scientific = FALSE)),
         median = ifelse(is.na(median), "ND",
                       format(round(median, digits = 1), scientific = FALSE))) %>% 
  pivot_wider(id_cols = CAS, 
              names_from = `Sample Date`,
              values_from = c("mean", "min", "max", "median", "n_dets", "samples")) %>% 
  arrange(match(CAS, df$CAS)) %>% 
  select(CAS,
         min_2010, max_2010, median_2010, mean_2010, n_dets_2010, samples_2010,
         min_2014, max_2014, median_2014, mean_2014, n_dets_2014, samples_2014)

writeData(wb, sheet = "Table S2",
          startRow = 5, startCol = 24 + ncol(df_tox), 
          x = select(chem_stats, -CAS))

saveWorkbook(wb, file = "test.xlsx", overwrite = TRUE)
saveWorkbook(wb, file = file.path(path_to_data, "Supplemental", "Table2_Combo.xlsx"), overwrite = TRUE)
