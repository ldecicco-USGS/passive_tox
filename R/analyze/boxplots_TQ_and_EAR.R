# Graph ECOTOX TQ values for all 

#### Setup ####
library(toxEval)
library(tidyverse)
#NOTE: Add path to path_to_file!!!
path_to_data <- Sys.getenv("PASSIVE_PATH")
path_to_file <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks.xlsx")
path_to_file_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx")
path_to_file_non_toxcast <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx")
path_to_file_EARs <- file.path(path_to_data, "data", "data_for_git_repo","clean", "passive.xlsx")

tox_list <- create_toxEval(path_to_file)
chemical_summary <- get_chemical_summary(tox_list)
benchmarks <- tox_list$benchmarks
######################################
bio_plot <- plot_tox_boxplots(chemical_summary, 
                              category = 'Chemical',
                              mean_logic = FALSE,
                              hit_threshold = c(0.01),
                              title = 'Max TQ per site (only detections)',
                              plot_ND = TRUE,
                              sum_logic = FALSE)
bio_plot

pdf("R/Analyze/Plots/TQ_boxplots_all_chems.pdf",height = 15)
bio_plot
dev.off()


p <- ggplot(chemical_summary,aes(x=EAR,y=chnm)) +
  geom_boxplot() + 
  scale_x_log10() + 
  facet_wrap(.~Bio_category)
       
p






# Only chems not in toxcast
tox_list <- create_toxEval(path_to_file_non_toxcast)
chemical_summary <- get_chemical_summary(tox_list)
######################################
bio_plot <- plot_tox_boxplots(chemical_summary, 
                              category = 'Chemical',
                              mean_logic = FALSE,
                              hit_threshold = c(0.01),
                              title = 'Max TQ per site (only detections)',
                              plot_ND = TRUE,
                              sum_logic = FALSE)
bio_plot

pdf("R/Analyze/Plots/TQ_boxplots_non_toxcast.pdf",height = 15)
bio_plot
dev.off()




# Only chems in toxcast
tox_list <- create_toxEval(path_to_file_toxcast)
chemical_summary <- get_chemical_summary(tox_list)
######################################
bio_plot <- plot_tox_boxplots(chemical_summary, 
                              category = 'Chemical',
                              mean_logic = FALSE,
                              hit_threshold = c(0.01),
                              title = 'Max TQ per site (only detections)',
                              plot_ND = TRUE,
                              sum_logic = FALSE)
bio_plot

pdf("R/Analyze/Plots/TQ_boxplots_toxcast.pdf",height = 15)
bio_plot
dev.off()



# Now with EAR values
tox_list <- create_toxEval(path_to_file_EARs)
chemical_summary <- get_chemical_summary(tox_list)
######################################
bio_plot <- plot_tox_boxplots(chemical_summary, 
                              category = 'Chemical',
                              mean_logic = FALSE,
                              hit_threshold = c(0.001),
                              title = 'Summing EARs of a sample, taking the max of each site',
                              plot_ND = TRUE)
bio_plot

pdf("R/Analyze/Plots/EAR_boxplots.pdf",height = 15)
bio_plot
dev.off()



# To save:
# Fiddle with height and width (in inches) for best results:
# Change file name extension to save as png.
# ggplot2::ggsave(bio_plot, file='boxplot.pdf',
#                        height = 9,
#                        width = 11)

tox_list <- create_toxEval(path_to_file)
chemical_summary <- get_chemical_summary(tox_list)
