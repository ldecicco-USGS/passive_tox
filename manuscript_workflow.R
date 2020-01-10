library(tidyverse)
library(toxEval)
library(openxlsx)

path_to_data <- Sys.getenv("PASSIVE_PATH")

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