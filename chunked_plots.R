library(toxEval)
library(ggplot2)
library(dplyr)

path_to_file <- 'cleanedData/passive.xlsx' 
tox_list <- create_toxEval(path_to_file)
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

chem_info <- tox_list$chem_info

plot_tox_boxplots(chemical_summary, "Chemical Class")

dir.create("plots", showWarnings = FALSE)

class_plot <- plot_tox_boxplots(chemical_summary, "Chemical Class")
ggsave(plot = class_plot, width = 11, height = 9,
       filename = "plots/class_boxplots.pdf")

chem_plot <- plot_tox_boxplots(chemical_summary, "Chemical")


ggsave(plot = chem_plot, width = 11, height = 33,
       filename = "plots/chem_boxplots.pdf")

pdf("plots/test.pdf", height = 9, width = 11)
for(i in levels(chemical_summary$Class)){
  
  sub_sum <- filter(chemical_summary, Class == i)
  sub_plot <- plot_tox_boxplots(sub_sum, "Chemical", title = i)
  print(sub_plot)
  # ggsave(plot = sub_plot, width = 11, height = 9,
  #        filename = paste0("plots/",i,"_boxplots.pdf"))
}
dev.off()

heat_plot <- plot_tox_heatmap(chemical_summary, chem_site = tox_list$chem_site, category = "Chemical")
ggsave(plot = heat_plot, width = 15, height = 33,
       filename = "plots/chem_heat.pdf")

pdf("plots/heat_chunked.pdf", height = 9, width = 11)
for(i in levels(chemical_summary$Class)){
  
  sub_sum <- filter(chemical_summary, Class == i)
  sub_sum$Class <- factor(sub_sum$Class, levels = i)
  sub_chems <- select(chemical_summary, chnm, Class) %>%
    distinct() %>%
    filter(Class == i) %>%
    pull(chnm)
  
  sub_sum$chnm <- factor(sub_sum$chnm, levels = as.character(sub_chems))
  
  sub_plot <- plot_tox_heatmap(sub_sum, title = i,
                               chem_site = tox_list$chem_site, category = "Chemical")
  print(sub_plot)
  # ggsave(plot = sub_plot, width = 11, height = 9,
  #        filename = paste0("plots/",i,"_boxplots.pdf"))
}
dev.off()
