library(toxEval)
library(ggplot2)

source("get_data.R")

chem_boxes <- plot_tox_boxplots(chemicalSummary, category = "Chemical")

dir.create("plots",showWarnings = FALSE)
ggsave(chem_boxes, filename = "plots/chemical_boxplots.pdf",height = 20,width = 9)

# Just pharms:
library(dplyr)
chem_sum_pharm <- filter(chemicalSummary, Class == "pharms")
pharm_boxes <- plot_tox_boxplots(chem_sum_pharm, category = "Chemical")
ggsave(pharm_boxes, filename = "plots/pharm_boxplots.pdf",height = 11,width = 9)

chem_sum_pharm_2 <- chem_sum_pharm 
chem_sum_pharm_2$Class <- as.character(chem_sum_pharm_2$date)
pharm_boxes <- plot_tox_boxplots(chem_sum_pharm_2, category = "Chemical")
ggsave(pharm_boxes, filename = "plots/pharm_boxplots_by_year.pdf",height = 11,width = 9)




