library(tidyverse)
library(toxEval)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(readxl)


create_triple_fig <- function(){
  tox_list <- create_toxEval(file.path(Sys.getenv("PASSIVE_PATH"),
                                       "data","data_for_git_repo","clean",
                                       "passive.xlsx"))
  
  source(file = "R/analyze/get_sites_ready.R")
  site_info <- prep_site_list(tox_list$chem_site)
  cas_final <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                 "data","data_for_git_repo","clean",
                                 "cas_df.rds"))
  
  tox_list_concentrations <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                            "data","data_for_git_repo","clean",
                                            "tox_list_concentrations.rds"))
  
  graphData_tox_det <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                         "data", "data_for_git_repo",
                                         "clean","graphData_tox_det.rds"))
  
  graphData_conc_det_match <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                                "data", "data_for_git_repo",
                                                "clean","graphData_conc_det_match.rds"))
   
  axis_num <- 6
  
  source(file = "R/report/combo_plot2.R")
  #########################################
  # # Get PCB?
  # source(file = "R/analyze/data_reader.R")
  # cas_df <- all_cas("data/raw/cas.xlsx")
  # x <- generic_file_opener(file_name = file.path(Sys.getenv("PASSIVE_PATH"),
  #                                                "data","data_for_git_repo","raw",
  #                                                "general_2014.xlsx"),
  #                          cas_df,
  #                     n_max = 45, 
  #                     sheet = "OC-PCB-PBDE",
  #                     site_sheet = "site info",
  #                     year = 2014)  %>% 
  #   filter(chnm %in% "Total Pcbs")
  #          
  # pcb <- data.frame(CAS = "1336-36-3",
  #                   chnm = "Total PCBs",
  #                   Class = "PCBs",
  #                   stringsAsFactors = FALSE)
  # 
  # tox_list_concentrations$chem_info <- bind_rows(tox_list_concentrations$chem_info,
  #                                                pcb)
  # 
  # data_pcb <- x %>% 
  #   select(SiteID, `Sample Date`, CAS, Value, comment )
  #   
  # tox_list_concentrations$chem_data <- bind_rows(tox_list_concentrations$chem_data,
  #                                                data_pcb)
  # tox_list_concentrations$benchmarks <- bind_rows(tox_list_concentrations$benchmarks,
  #                                                 mutate(pcb,
  #                                                        endPoint = "Concentration",
  #                                                        ACC_value = 1,
  #                                                        groupCol = "Concentration"))
  # cas_final <- bind_rows(cas_final,
  #                        pcb)
  
  chemicalSummary_conc <- get_chemical_summary(tox_list_concentrations)
  
  #######################
  chemicalSummary_conc_no_match = chemicalSummary_conc %>%
    filter(!(CAS %in% unique(graphData_tox_det$CAS)),
           EAR > 0)
  
  graphData_conc_no_match = graph_chem_data_CAS(chemicalSummary_conc_no_match) %>%
    mutate(guide_side = "Concentration [\U003BCg/L]") %>%
    left_join(select(cas_final, CAS, chnm), by="CAS")
  
  full_classes <- c(levels(graphData_tox_det$Class),
                    levels(graphData_conc_no_match$Class)[!(levels(graphData_conc_no_match$Class) %in% levels(graphData_tox_det$Class))])
  
  graphData_tox_det$Class <- factor(as.character(graphData_tox_det$Class), levels = full_classes)
  graphData_conc_no_match$Class <- factor(as.character(graphData_conc_no_match$Class), levels = full_classes)
  
  matches <- fancy_combo(graphData_tox_det, 
                         graphData_conc_det_match, 
                         tox_list, 
                         axis_size = axis_num)
  
  n_chems_matches <- length(unique(graphData_tox_det$chnm))
  
  graphData_empty <- graphData_conc_no_match[FALSE,]
  
  gd_no_match <- combine_gd(graphData_conc_no_match, graphData_empty)
  
  n_chems_no_match <- length(unique(gd_no_match$chnm))
  
  color_map <- class_colors(chemicalSummary_conc)
  toxPlot_no_match <- combo_plot_matches_2(gd_no_match,
                                     axis_size = axis_num,
                                     color_map)
  text_df_c <- label_info(gd_no_match, labels_to_use = "C")
  toxPlot_no_match_w_lab <- add_label(toxPlot_no_match, text_df_c)
  no_axis_no_match <- strip_graph(toxPlot_no_match_w_lab)
  site_counts_df_no_match <- site_counts(tox_list_concentrations$chem_data, no_axis_no_match$data)
  site_graph_no_match <- site_count_plot(site_counts_df_no_match,
                                         axis_size = axis_num)
  
  l2 <- get_legend(toxPlot_no_match)

  plot_out <- plot_grid(
    matches$site_graph,
    matches$no_axis,
    plot_grid(
      plot_grid(
        site_graph_no_match, 
        no_axis_no_match,
        ncol = 2,
        rel_widths = c(2.25,3)
      ),
      plot_grid(
        l2,
        NULL,
        ncol=1
      ),
      nrow = 2, ncol = 1,
      rel_heights = c(n_chems_no_match,n_chems_matches-n_chems_no_match)
    ),
    rel_widths = c(2.75,4,5),
    nrow=1,ncol=3
  )
  
  return(plot_out)
  
}

# pdf(file.path(Sys.getenv("PASSIVE_PATH"),"Figures/Polished figures/triple_graph.pdf"),
#               width = 9, height = 11, onefile=FALSE)
# create_triple_fig()
# dev.off()
