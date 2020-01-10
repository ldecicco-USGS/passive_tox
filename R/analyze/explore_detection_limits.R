library(toxEval)
library(tidyverse)

create_DL_plot <- function(tox_list){
  
  chem_data_det <- tox_list$chem_info %>% 
    select(-chnm, -Class, -sites_det, -sites_tested) %>% 
    pivot_longer(cols = c(-CAS), 
                 names_to = "DL",
                 values_to = "Value") %>% 
    filter(!is.na(Value)) %>% 
    mutate(`Sample Date` = 1,
           SiteID = tox_list$chem_site$SiteID[1])
  
  tox_list_dl <- tox_list
  tox_list_dl$chem_data <- chem_data_det
  
  ACC <- get_ACC(tox_list_dl$chem_info$CAS)
  ACC <- remove_flags(ACC)
  
  cleaned_ep <- clean_endPoint_info(end_point_info)
  filtered_ep <- filter_groups(cleaned_ep)
  
  chemical_summary <- get_chemical_summary(tox_list_dl, ACC, filtered_ep)
  
  ch_plot <- plot_tox_boxplots(chemical_summary, 
                               "Chemical", x_label = "EARs of detection limits")
  return(ch_plot)

}  
# ggsave(ch_plot, filename = "plots/det_limits.pdf", height = 22, width = 11)
