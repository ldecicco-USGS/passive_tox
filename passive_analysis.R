library(drake)
library(tidyverse)
library(toxEval)
library(cowplot)

options(drake_make_menu = FALSE)

source(file = "passive_data_setup.R")
source(file = "R/report/stack_plots.R")
source(file = "R/report/combo_graph_function.R")
source(file = "R/analyze/graph_chem_data_CAS.R")
source(file = "R/report/plot_tox_endpoints_manuscript.R")

loadd(cas_df)

cas_final =  cas_df %>%
  filter(!duplicated(CAS)) %>%
  mutate(chnm = tools::toTitleCase(chnm))

data_analysis_plan <- drake_plan(

  tox_list = create_toxEval(file_in("data/clean/passive.xlsx")),
  ACC = get_ACC(tox_list$chem_info$CAS) %>%
    remove_flags(),
  cleaned_ep = clean_endPoint_info(end_point_info),
  filtered_ep = filter_groups(cleaned_ep, 
                               groupCol = 'intended_target_family',
                               assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                               remove_groups = c('Background Measurement','Undefined')),
  
  chemicalSummary = get_chemical_summary(tox_list, ACC, filtered_ep),
  benchmarks = tox_list$chem_data %>%
    select(CAS) %>%
    distinct() %>%
    left_join(cas_df, by="CAS") %>%
    mutate(endPoint = "Concentration",
           Value = 1,
           groupCol = "Concentration") %>%
    data.frame(),
  tox_list_concentrations = as.toxEval(tox_list, benchmarks = benchmarks),
  chemicalSummary_conc = get_chemical_summary(tox_list_concentrations) %>%
    distinct() %>%
    filter(!is.na(CAS)),
  # Think about duplicate cas's here!
  graphData_tox = graph_chem_data_CAS(chemicalSummary) %>%
    mutate(guide_side = "ToxCast",
           guide_up = "A"),
  graphData_conc = graph_chem_data_CAS(chemicalSummary_conc) %>%
    mutate(guide_side = "Concentration",
           guide_up = "A"),
  toxPlot_ear_conc = combo_plot_matches(gd_1 = graphData_tox,
                                         gd_2 = graphData_conc,
                                        cas_key = cas_final,
                                         thres_1 = NA,
                                         thres_2 = NA,
                                         drop = FALSE),
  ggsave(toxPlot_ear_conc, 
         filename = file_out("plots/EAR_Conc_all.pdf"),width = 9, height = 22),
  graphData_conc_det = filter(graphData_conc, meanEAR > 0),
  graphData_tox_det = filter(graphData_tox, meanEAR > 0),
  toxPlot_ear_conc_detects = combo_plot_matches(gd_1 = graphData_tox_det,
                                                gd_2 = graphData_conc_det,
                                                thres_1 = NA,
                                                thres_2 = NA,
                                                drop = FALSE,
                                                cas_key = cas_final),
  ggsave(toxPlot_ear_conc_detects, 
         filename = file_out("plots/EAR_Conc_detects.pdf"),width = 9, height = 22),
  graphData_conc_det_match = filter(graphData_conc_det, CAS %in% unique(graphData_tox_det$CAS)),
  toxPlot_ear_conc_matches = combo_plot_matches(gd_1 = graphData_tox_det,
                                        gd_2 = graphData_conc_det_match,
                                        thres_1 = NA,
                                        thres_2 = NA,
                                        drop = FALSE,
                                        cas_key = cas_final),
  ggsave(toxPlot_ear_conc_matches, 
         filename = file_out("plots/EAR_Conc_only_detects_matches.pdf"),width = 9, height = 22),
  AOP = readr::read_csv(file_in(readd(AOP_crosswalk))) %>%
                          select(endPoint=`Component Endpoint Name`, ID=`AOP #`) %>%
                          distinct(),
  aop_graph = plot_tox_endpoints_manuscript(chemicalSummary, AOP,
                                            category = "Chemical", 
                                            font_size = 7,title = " ",
                                            pallette = c("steelblue", "white")),
  save_plot(filename = file_out("plots/AOPs.pdf"),
            plot = aop_graph, base_width = 11),
  save_plot(filename = file_out("plots/AOPs.png"),
            plot = aop_graph, base_width = 11),
  
  site_info = prep_site_list(tox_list$chem_site),
  stack_plot = plot_tox_stacks_manuscript(chemicalSummary, site_info, 
                                           category = "Chemical Class"),
  save_plot(filename = file_out("plots/stacks.png"),
            plot = stack_plot, base_width = 11, base_height = 7),
  save_plot(filename = file_out("plots/stacks.pdf"),
            plot = stack_plot, base_width = 11, base_height = 7)

)


drake_config(data_analysis_plan)
# config <- drake_config(data_analysis_plan)
# vis_drake_graph(config, build_times = "none")

