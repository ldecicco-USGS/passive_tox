library(tidyverse)
library(toxEval)
library(openxlsx)

path_to_data <- Sys.getenv("PASSIVE_PATH")

# shiny::runApp("apps/Mixture_Exploration/")

# Data setup is done in drake plan
# It's better to run separately, and shouldn't
# be needed to run this script

source(file = "read_chemicalSummary.R")

unique(tox_list$chem_data$CAS)[which(!(unique(tox_list$chem_data$CAS)) %in% tox_list$chem_info$CAS)]

tox_list$exclusions <- tox_list$exclusions %>% 
  filter(!is.na(CAS) & !is.na(endPoint))
###################################################
# Mixtures stuff:
EAR_thresh <- 0.001
site_thresh_percent <- 10

source(file = "R/mixtures/mix_script.R")
source(file = "R/mixtures/prepare_mixture_data.R")
site_thresh <- ceiling((site_thresh_percent/100) * nrow(tox_list$chem_site))

mix_df <- get_final_mixtures(chemicalSummary,
                             EAR_thresh,
                             site_thresh)
dir.create("R/mixtures/out", showWarnings = FALSE)
saveRDS(mix_df,"R/mixtures/out/mixtures_table.rds")


wb <- create_Excel_wb_mix(mix_df)

saveWorkbook(wb, file.path(Sys.getenv("PASSIVE_PATH"),"Tables",
                           "Mixtures.xlsx"), overwrite = TRUE)

################################################

# Figure 1 = Map
# Figure 2 = triple_graph.pdf
# Figure 3 = new_stack_w_table.pdf
rmarkdown::render("manuscript_figures.Rmd", 
                  output_dir = file.path(Sys.getenv("PASSIVE_PATH"),
                                         "Figures"))


# Figure SI A-D = 
# Run the stacked_supplement_w_captions.Rmd file
rmarkdown::render("stacked_supplement_w_captions.Rmd", 
       output_dir = file.path(Sys.getenv("PASSIVE_PATH"),
                              "Supplemental"))
# This graph includes chemicals measured below detection levels...
# which is why it is longer and there appear to be more endpoints
# than the numbers we wrote in the text (those were all for detected chemicals)

# Create SI tables:

df_assays <- chemicalSummary %>% 
  select(endPoint) %>% 
  distinct() %>% 
  left_join(select(end_point_info, 
                   endPoint = assay_component_endpoint_name,
                   source = assay_source_long_name), by="endPoint") %>% 
  rename(`ToxCast Assay` = endPoint,
         `Assay Source` = source)

assays_left <- chemicalSummary %>% 
  select(CAS, endPoint) %>% 
  distinct() %>% 
  group_by(CAS) %>% 
  summarize(`Assays in study` = length(unique(endPoint)))

AOP_crosswalk <- data.table::fread(file.path(path_to_data, "data/data_for_git_repo/raw", "AOP_crosswalk.csv"), data.table = FALSE)

AOP_sup <- AOP_crosswalk %>% 
  select(AOP_shortname, `AOP #`, `AOP Title`, `KE#`, `Key Event Name`) %>% 
  distinct() %>% 
  filter(AOP_shortname != "")

all_join <- join_everything_fnx(chemicalSummary) %>% 
  separate(genes, sep  = ",", into = letters[1:5]) %>% 
  pivot_longer(letters[1:5], "names", values_to = "Gene") %>% 
  select(-names) %>% 
  filter(!is.na(Gene),
         Gene != "") %>% 
  mutate(AOPs = ifelse(AOPs == "", "No", "Yes"),
         Panther = ifelse(pathways == "", "No", "Yes")) %>% 
  select(-pathways)

AOP_pan <- chemicalSummary %>% 
  select(chnm, endPoint) %>% 
  distinct() %>% 
  left_join(all_join, by = "endPoint") %>% 
  select(Gene, chnm, endPoint, AOPs, Panther) %>% 
  arrange(Gene)

#Run ECOTOX analysis
suppressMessages(source("R/Analyze/ECOTOX_workflow.R"))

#Combine EAR and TQ-based chem priorities
source("R/Analyze/combine_TQ_EAR_chem_priorities.R")

################################################
#Create table 2: summary of priority chemicals
source("R/report/chem_priority_summary_Table2.R")

table_2 <- get_table_2()
names(table_2) <- c("Chemical use class","Chemical name","Sites monitored","Individual","Mixtures"," ","Group 1","Group 2")

ms_tables_wb <- loadWorkbook(file.path(path_to_data,"Tables","Manuscript_tables.xlsx"))
writeData(wb = ms_tables_wb,sheet = "Table 2. Chemical Summary",x = table_2,startCol = 1,startRow = 4)
saveWorkbook(wb = ms_tables_wb,file = file.path(path_to_data,"Tables","Manuscript_tables.xlsx"),overwrite = TRUE)

################################################
# Create the supplemental:
source(file = "R/analyze/update_SI_2.R")

tab_names <- c("SI-1 Site Table",
               "SI-2 Chemical Table",
               "SI-3 ToxCast Assays",
               "SI-4 ToxCast Exclusions",
               "SI-5 ECOTOX data",
               "SI-6 ECOTOX summary",
               "SI-7 AOP",
               "SI-8 AOP and Panther"
)

captions <- c("Table SI-1: Sampling locations and primary land cover characteristics for Great Lakes tributaries sampled by passive samplers, 2010-2014.",
              "Table SI-2: Chemicals, prescribed chemical use, passive sampler type, sampling rate, ToxCast assay information, detection levels, and summary statistics from resulting data for Great Lakes tributaries, 2010-2014",
              "Table SI-3: Toxcast assays included in data analysis for organic chemicals detected in Great Lakes tributary passive samplers, 2010-2014",
              "Table SI-4: ToxCast assays excluded from analysis due low quality dose-response curves based on anomalous values or lack of response ",
              "Table SI-5: Descriptive information for bioassay endpoints from the ECOTOX knowledgebase included in data analysis for organic chemicals detected in Great Lakes tributary passive samplers, 2010-2014", 
              "Table SI-6: Bioassay endpoint information for individual chemicals from the ECOTOX knowledgebase included in data analysis for organic chemicals detected in Great Lakes tributary passive samplers, 2010-2014. [\U003BCg/L: Units for min endpiont, max endpoint, and 5th centile of endpoints",
              "Table SI-7: Adverse outcome pathway (AOP) identification and key event information from the AOP Knowledgebase for AOPs relevant to  organic chemicals detected in Great Lakes tributary passive samplers, 2010-2014",
              "Table SI-8: ToxCast endpoints, gene targets, and association with adverse outcome pathways and gene ontology pathways (Panther) for individual chemicals from the ECOTOX knowledgebase included in data analysis for organic chemicals detected in Great Lakes tributary passive samplers, 2010-2014")

wb <- update_table_2(path_to_data, tab_names[2])

# SI-1: Site Table
addWorksheet(wb, tab_names[1])
header_st <- createStyle(textDecoration = "Bold")
writeData(wb = wb, sheet =  tab_names[1], colNames = FALSE, rowNames = FALSE,
          x = captions[1])
writeData(wb = wb, sheet = tab_names[1], startRow = 3,
          x = tox_list$chem_site, headerStyle = header_st)

worksheetOrder(wb) <- c(2,1)
#SI-2: Chemical Table
# See update_table_2


# SI-3: ToxCast Assays
addWorksheet(wb, tab_names[3])
writeData(wb = wb, sheet = tab_names[3], colNames = FALSE, rowNames = FALSE,
          x = captions[3])
writeData(wb = wb, sheet = tab_names[3], startRow = 3,
          x = df_assays, headerStyle = header_st)

# SI-4: Exclusions
addWorksheet(wb, tab_names[4])
writeData(wb = wb, sheet = tab_names[4], colNames = FALSE, rowNames = FALSE,
          x = captions[4])
writeData(wb = wb, sheet = tab_names[4], startRow = 3,
          x = rename(select(tox_list$exclusions, -Chemical),
                     `ToxCast Assay` = endPoint, Chemical = chnm), headerStyle = header_st)
#SI-5 ECOTOX data
addWorksheet(wb, tab_names[5])
writeData(wb = wb, sheet = tab_names[5], colNames = FALSE, rowNames = FALSE,
          x = captions[5])
writeData(wb = wb, sheet = tab_names[5], startRow = 3,
          x = tox_fw)

#SI-6 ECOTOX summary stats
addWorksheet(wb, tab_names[6])
writeData(wb = wb, sheet = tab_names[6], colNames = FALSE, rowNames = FALSE,
          x = captions[6])
writeData(wb = wb, sheet = tab_names[6], startRow = 3,
          x = tox_stats)

#SI-7 AOP crosswalk
addWorksheet(wb, tab_names[7])
writeData(wb = wb, sheet = tab_names[7], colNames = FALSE, rowNames = FALSE,
          x = captions[7])
writeData(wb = wb, sheet = tab_names[7], startRow = 3,
          x = AOP_sup,
          headerStyle = header_st)

#SI-8 AOP and Panther
addWorksheet(wb, tab_names[8])
writeData(wb = wb, sheet = tab_names[8], colNames = FALSE, rowNames = FALSE,
          x = captions[8])
writeData(wb = wb, sheet = tab_names[8], startRow = 3,
          x = AOP_pan,
          headerStyle = header_st)



# Save the whole thing:
saveWorkbook(wb, file.path(Sys.getenv("PASSIVE_PATH"),"Supplemental",
                           "Supplemental_32.xlsx"), overwrite = TRUE)

