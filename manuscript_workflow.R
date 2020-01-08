library(tidyverse)
library(toxEval)

path_to_data <- Sys.getenv("PASSIVE_PATH")

# Data setup is done in drake plan
# It's better to run separately, and shouldn't
# be needed to run this script

source(file = "read_chemicalSummary.R")
tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                     "data","data_for_git_repo","clean",
                                     "passive.xlsx"))

site_info <- tox_list$chem_site

EAR_thresh <- 0.001
site_thresh_percent <- 10
site_thresh <- ceiling((site_thresh_percent/100) * nrow(site_info))

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
