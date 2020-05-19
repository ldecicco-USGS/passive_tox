library(drake)
library(tidyverse)
library(toxEval)
library(ggpubr)

options(drake_make_menu = FALSE)

path_to_data <- Sys.getenv("PASSIVE_PATH")

# source(file = "passive_data_setup.R")
source(file = "R/report/stack_plots.R")
source(file = "R/report/combo_plot2.R")
source(file = "R/analyze/graph_chem_data_CAS.R")
source(file = "R/report/plot_tox_endpoints_manuscript.R")

data_analysis_plan <- drake_plan(
  cas_final = readRDS(file_in(!!file.path(path_to_data, "data/data_for_git_repo/clean/cas_df.rds"))),
  tox_list = create_toxEval(file_in(!!file.path(path_to_data, "data/data_for_git_repo/clean/passive.xlsx"))),
  ACC = get_ACC(tox_list$chem_info$CAS) %>%
    remove_flags(),
  cleaned_ep = clean_endPoint_info(end_point_info),
  filtered_ep = filter_groups(cleaned_ep, 
                              groupCol = 'intended_target_family',
                              remove_groups = c('Background Measurement','Undefined','Cell Cycle','NA')),
  
  chemicalSummary = get_chemical_summary(tox_list, ACC, filtered_ep),
  chem_sum_save_sync = saveRDS(chemicalSummary,
                          file = file_out(!!file.path(path_to_data,"data/data_for_git_repo/clean/chemical_summary.rds"))),
  benchmarks = tox_list$chem_data %>%
    select(CAS) %>%
    distinct() %>%
    left_join(select(cas_final, CAS, Class, chnm), by="CAS") %>%
    mutate(endPoint = "Concentration",
           Value = 1,
           groupCol = "Concentration") %>%
    data.frame(),
  tox_list_concentrations = as.toxEval(tox_list, benchmarks = benchmarks),
  tox_list_concentrations_sync = saveRDS(tox_list_concentrations, 
                                         file = file_out(!!file.path(path_to_data,
                                                                     "data/data_for_git_repo/clean/tox_list_concentrations.rds"))),
  chemicalSummary_conc = get_chemical_summary(tox_list_concentrations) %>%
    distinct() %>%
    filter(!is.na(CAS)),
  chemicalSummary_conc_save = saveRDS(chemicalSummary_conc, 
                                      file = file_out(!!file.path(path_to_data,"data/data_for_git_repo/clean/chemicalSummary_conc.rds"))),
  graphData_tox = graph_chem_data_CAS(chemicalSummary) %>%
    mutate(guide_side = "ToxCast [EAR]") %>%
    left_join(select(cas_final, CAS, chnm), by="CAS"),
  graphData_conc = graph_chem_data_CAS(chemicalSummary_conc) %>%
    mutate(guide_side = "Concentration [\U003BCg/L]") %>%
    left_join(select(cas_final, CAS, chnm), by="CAS"),
  toxPlot_ear_conc = fancy_combo(graphData_tox,
                                 graphData_conc,
                                 tox_list),
  chemicalSummary_conc_det = filter(chemicalSummary_conc, EAR > 0),
  chemicalSummary_tox_det = filter(chemicalSummary, EAR > 0),
  graphData_conc_det = graph_chem_data_CAS(chemicalSummary_conc_det) %>%
    mutate(guide_side = "Concentration [\U003BCg/L]") %>%
    left_join(select(cas_final, CAS, chnm), by="CAS"),
  graphData_tox_det = graph_chem_data_CAS(chemicalSummary_tox_det) %>%
    mutate(guide_side = "ToxCast [EAR]") %>%
    left_join(select(cas_final, CAS, chnm), by="CAS"),
  graphData_tox_det_out = saveRDS(graphData_tox_det, 
                                  file = file_out(!!file.path(path_to_data,
                                                              "data/data_for_git_repo/clean/graphData_tox_det.rds"))),
  toxPlot_ear_conc_detects = fancy_combo(graphData_tox_det,
                                         graphData_conc_det,
                                         tox_list),
  chemicalSummary_conc_det_match = filter(chemicalSummary_conc_det, CAS %in% unique(graphData_tox_det$CAS)),
  graphData_conc_det_match = graph_chem_data_CAS(chemicalSummary_conc_det_match) %>%
    mutate(guide_side = "Concentration [\U003BCg/L]") %>%
    left_join(select(cas_final, CAS, chnm), by="CAS"),
  graphData_conc_det_match_save = saveRDS(graphData_conc_det_match, 
                                          file = file_out(!!file.path(path_to_data,
                                                                      "data/data_for_git_repo/clean/graphData_conc_det_match.rds"))),
  toxPlot_ear_conc_matches = fancy_combo(graphData_tox_det,
                                         graphData_conc_det_match,
                                         tox_list),
  graphData_conc_det_match_filter = graphData_conc_det_match,
  graphData_tox_det_filter = graphData_tox_det,
  toxPlot_ear_conc_matches_filter = fancy_combo(graphData_tox_det_filter,
                                                graphData_conc_det_match_filter,
                                                tox_list),
  AOP = readr::read_csv(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/AOP_crosswalk.csv"))) %>%
                          select(endPoint=`Component Endpoint Name`, ID=`AOP #`) %>%
                          distinct(),
  aop_graph = plot_tox_endpoints_manuscript(chemicalSummary, AOP,
                                            category = "Chemical", 
                                            font_size = 7,title = " ",
                                            pallette = c("steelblue", "white")),
  site_info = prep_site_list(tox_list$chem_site), #this makes the data frame factors with the right order
  stack_plot = plot_tox_stacks_manuscript(chemicalSummary, 
                                          site_info, 
                                          font_size = 6,
                                          category = "Chemical Class")

)

vis_drake_graph(data_analysis_plan, build_times = "none")

make(data_analysis_plan, trigger = trigger(condition=TRUE))

# drake_config(data_analysis_plan)
# In R console:
# r_make("passive_data_setup.R")


# loadd(aop_graph)
# pdf("plots/AOP_v3.pdf", width = 4.5, height = 4.5)
# ggarrange(
#   aop_graph$count_plot, aop_graph$stackedPlot,
#   aop_graph$chem_plot, aop_graph$aop_plot,nrow = 1,ncol = 4,
#   widths =  c(4/10,4/10,1/10,1/10)
# )
# dev.off()
# 
# loadd(toxPlot_ear_conc)
# pdf("plots/EAR_Conc_all.pdf", width = 4.5, height = 22, onefile=FALSE)
# ggarrange(
#   toxPlot_ear_conc$site_graph, 
#   toxPlot_ear_conc$no_axis,
#   widths =  c(3.25/9, 5.75/9),
#   common.legend = TRUE, legend = "bottom"
# )
# dev.off()
# 
# loadd(toxPlot_ear_conc_detects)
# pdf("plots/EAR_Conc_detects.pdf", width = 4.5, height = 11, onefile=FALSE)
# ggarrange(
#   toxPlot_ear_conc_detects$site_graph, 
#   toxPlot_ear_conc_detects$no_axis,
#   widths =  c(3.25/9, 5.75/9),
#   common.legend = TRUE, legend = "bottom"
# )
# dev.off()
# 
# loadd(toxPlot_ear_conc_matches)
# pdf("plots/EAR_Conc_detects_match.pdf", width = 4.75, height = 9, onefile=FALSE)
# ggarrange(
#   toxPlot_ear_conc_matches$site_graph, 
#   toxPlot_ear_conc_matches$no_axis,
#   widths =  c(3.25/9, 5.75/9),
#   common.legend = TRUE, legend = "bottom"
# )
# dev.off()
# 
# 
# ###########################
# loadd(graphData_tox_det)
# loadd(graphData_conc_det_match)
# loadd(tox_list)
# source(file = "R/report/combo_plot2.R")
# toxPlot_ear_conc_matches = fancy_combo(graphData_tox_det,
#                                        graphData_conc_det_match,
#                                        tox_list)
# 
# pdf("plots/EAR_Conc_detects_match.pdf", width = 4.5, height = 9, onefile=FALSE)
# ggarrange(
#   toxPlot_ear_conc_matches$site_graph, 
#   toxPlot_ear_conc_matches$no_axis,
#   widths =  c(3.25/9, 5.75/9),
#   common.legend = TRUE, legend = "bottom"
# )
# dev.off()
