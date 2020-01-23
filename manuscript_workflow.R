library(tidyverse)
library(toxEval)
library(openxlsx)

path_to_data <- Sys.getenv("PASSIVE_PATH")

# shiny::runApp("apps/Mixture_Exploration/")

# Data setup is done in drake plan
# It's better to run separately, and shouldn't
# be needed to run this script

source(file = "read_chemicalSummary.R")
tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))


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


ALL_TOX_DATA <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                  "data","data_for_git_repo","raw",
                                  "all_tox.rds"))

chem_data <- tox_list$chem_info

ALL_TOX_DATA_in_study <-  ALL_TOX_DATA %>% 
  select(CAS = casn, endPoint=aenm, modl_acc, flags, hitc) %>% 
  filter(CAS %in% chem_data$CAS) %>% 
  group_by(CAS) %>% 
  summarize(`Total ToxCast assays` = length(unique(endPoint)),
            `Assays with hits` = length(unique(endPoint[hitc == 1])))

assays_left <- chemicalSummary %>% 
  select(CAS, endPoint) %>% 
  distinct() %>% 
  group_by(CAS) %>% 
  summarize(`Assays in study` = length(unique(endPoint)))

chem_data <- chem_data %>% 
  left_join(ALL_TOX_DATA_in_study, by="CAS") %>% 
  left_join(assays_left, by="CAS")

AOP_crosswalk <- data.table::fread(file.path(path_to_data, "data/data_for_git_repo/raw", "AOP_crosswalk.csv"), data.table = FALSE)

################################################
# Create the supplemental:
wb <- createWorkbook()
# SI-1: Site Table
addWorksheet(wb, "SI-1 Site Table")
header_st <- createStyle(textDecoration = "Bold")
writeData(wb = wb, sheet =  "SI-1 Site Table", colNames = FALSE, rowNames = FALSE,
          x = "Table SI-1: Site information")
writeData(wb = wb, sheet =  "SI-1 Site Table", startRow = 3,
          x = tox_list$chem_site, headerStyle = header_st)

#SI-2: Chemical Table
addWorksheet(wb, "SI-2 Chemical Table")
writeData(wb = wb, sheet =  "SI-2 Chemical Table", colNames = FALSE, rowNames = FALSE,
          x = "Table SI-2: Chemical Information")
writeData(wb = wb, sheet =  "SI-2 Chemical Table", startRow = 3,
          x = chem_data, headerStyle = header_st)

# SI-3: POCIS sampling rates
addWorksheet(wb, "SI-3 POCIS sampling rates")
writeData(wb = wb, sheet =  "SI-3 POCIS sampling rates", colNames = FALSE, rowNames = FALSE,
          x = "Table SI-3: POCIS sampling rates")


# SI-4: ToxCast Assays
addWorksheet(wb, "SI-4 ToxCast Assays")
writeData(wb = wb, sheet =  "SI-4 ToxCast Assays", colNames = FALSE, rowNames = FALSE,
          x = "Table SI-4: ToxCast Assays")
writeData(wb = wb, sheet =  "SI-4 ToxCast Assays", startRow = 3,
          x = df_assays, headerStyle = header_st)

# SI-5: Exclusions
addWorksheet(wb, "SI-5 Exclusions")
writeData(wb = wb, sheet =  "SI-5 Exclusions", colNames = FALSE, rowNames = FALSE,
          x = "Table SI-5: Exclusions")
writeData(wb = wb, sheet =  "SI-5 Exclusions", startRow = 3,
          x = rename(tox_list$exclusions,
                     `ToxCast Assay` = endPoint, Chemical = chnm), headerStyle = header_st)

#SI-6 AOP crosswalk
addWorksheet(wb, "SI-6 AOP")
writeData(wb = wb, sheet =  "SI-6 AOP", colNames = FALSE, rowNames = FALSE,
          x = "Table SI-6: AOP")
writeData(wb = wb, sheet =   "SI-6 AOP", startRow = 3,
          x = mtcars,
          headerStyle = header_st)

# Save the whole thing:
saveWorkbook(wb, file.path(Sys.getenv("PASSIVE_PATH"),"Supplemental",
                           "Supplemental.xlsx"), overwrite = TRUE)

