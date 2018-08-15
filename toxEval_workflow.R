library(toxEval)
library(ggplot2)

tox_list <- create_toxEval("cleanedData/passive.xlsx")

ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info = end_point_info)
filtered_ep <- filter_groups(cleaned_ep)

lakes_ordered <- c("Lake Superior",
                   "Lake Michigan",
                   "Lake Huron",
                   "Lake Erie",
                   "Lake Ontario")

tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping,
                                           levels=lakes_ordered)

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

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
