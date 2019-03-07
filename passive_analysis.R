library(drake)
library(tidyverse)
library(toxEval)

source("R/report/combo_graph_function.R")
source("passive_data_setup.R")
loadd(cas_df)

data_analysis_plan <- drake_plan(
  tox_list = create_toxEval(file_out_data),
  
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
  chemicalSummary_conc = get_chemical_summary(tox_list_concentrations),
  # Think about duplicate cas's here!
  graphData_tox = graph_chem_data(chemicalSummary) %>%
    mutate(guide_side = "ToxCast",
           guide_up = "A") %>%
    left_join(distinct(select(chemicalSummary, CAS, chnm)), by="chnm") %>%
    select(-chnm) %>%
    left_join(cas_df, by = "CAS"),
  graphData_conc = graph_chem_data(chemicalSummary_conc) %>%
    mutate(guide_side = "Concentration",
           guide_up = "A"),
  toxPlot_ear_conc = combo_plot_matches(gd_1 = graphData_tox,
                                         gd_2 = graphData_conc,
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
                                                drop = FALSE),
  ggsave(toxPlot_ear_conc_detects, 
         filename = file_out("plots/EAR_Conc_detects.pdf"),width = 9, height = 22),
  toxPlot_ear_conc_matches = combo_plot_matches(gd_1 = graphData_tox_det,
                                        gd_2 = graphData_conc_det,
                                        thres_1 = NA,
                                        thres_2 = NA,
                                        drop = TRUE),
  ggsave(toxPlot_ear_conc_matches, 
         filename = file_out("plots/EAR_Conc_only_detects_matches.pdf"),width = 9, height = 22)
  
  
)

config <- drake_config(data_analysis_plan)
vis_drake_graph(config, build_times = "none")
make(data_analysis_plan)
loadd(toxPlot_ear_conc)
toxPlot_ear_conc
