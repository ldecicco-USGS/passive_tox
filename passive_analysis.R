library(drake)
library(tidyverse)
library(toxEval)
library(ggpubr)

options(drake_make_menu = FALSE)

source(file = "passive_data_setup.R")
source(file = "R/report/stack_plots.R")
source(file = "R/report/combo_plot2.R")
source(file = "R/analyze/graph_chem_data_CAS.R")
source(file = "R/report/plot_tox_endpoints_manuscript.R")

loadd(cas_df)

cas_final =  cas_df %>%
  filter(!duplicated(CAS)) %>%
  mutate(chnm = tools::toTitleCase(chnm))

cas_final$chnm[cas_final$chnm == "Deet"] <- "DEET"
cas_final$chnm[cas_final$chnm == "O,p'-Ddd"] <- "o,p'-DDD"
cas_final$chnm[cas_final$chnm == "P,p'-Ddd"] <- "p,p'-DDD"
cas_final$chnm[cas_final$chnm == "Pentachloroanisole (Pca)"] <- "PCA"
cas_final$chnm[cas_final$chnm == "Tributyl Phosphate (Tbp)"] <- "TBP"
cas_final$chnm[cas_final$chnm == "Hydrochlorothiazide (Hctz)"] <- "HCTZ"
cas_final$chnm[cas_final$chnm == "Tris(2−Chloroethyl)Phosphate (Tcep)"] <- "TCEP"
cas_final$chnm[cas_final$chnm == "O,p'−Ddt"] <- "o,p'−DDT"
cas_final$chnm[cas_final$chnm == "P,p'−Dde"] <- "p,p'−DDE"
cas_final$chnm[cas_final$chnm == "P,p'−Ddt"] <- "p,p'−DDT"
cas_final$chnm[cas_final$chnm == "O,p'−Dde"] <- "o,p'−DDE"
cas_final$chnm[cas_final$chnm == "Tris(1-Chloro-2-Propyl)Phosphate (Tcpp)"] <- "TCPP"
cas_final$chnm[cas_final$chnm == "Hexachlorobenzene (Hcb)"] <- "HCB"
cas_final$chnm[cas_final$CAS == "77-93-0"] <- "Triethyl Citrate "
cas_final$chnm[cas_final$CAS == "30306-93-5"] <- "Ethyl Citrate"
cas_final$chnm[grep("Pbde-", cas_final$chnm)] <- gsub(pattern = "Pbde-",
                                                      replacement = "PBDE-",
                                                      cas_final$chnm[grep("Pbde-", cas_final$chnm)])
cas_final <- rbind(cas_final, data.frame(CAS="34841-39-9",
                                         chnm="Bupropion",
                                         stringsAsFactors = FALSE))
cas_final$chnm[cas_final$CAS == "34911-55-2"] <- "Bupropion hydrochloride"

cas_final$chnm[grep(pattern = "Delta-Benzenehexachloride",cas_final$chnm)] <- "delta-Bhc"
cas_final$chnm[grep(pattern = "Beta-Benzenehexachloride",cas_final$chnm)] <- "beta-Bhc"
cas_final$chnm[grep(pattern = "Alpha-Benzenehexachloride", cas_final$chnm)] <- "alpha-Bhc"

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
  graphData_tox = graph_chem_data_CAS(chemicalSummary) %>%
    mutate(guide_side = "ToxCast [EAR]") %>%
    left_join(cas_final, by="CAS"),
  graphData_conc = graph_chem_data_CAS(chemicalSummary_conc) %>%
    mutate(guide_side = "Concentration [\U003BCg/L]") %>%
    left_join(cas_final, by="CAS"),
  toxPlot_ear_conc = fancy_combo(graphData_tox,
                                 graphData_conc,
                                 tox_list),
  chemicalSummary_conc_det = filter(chemicalSummary_conc, EAR > 0),
  chemicalSummary_tox_det = filter(chemicalSummary, EAR > 0),
  graphData_conc_det = graph_chem_data_CAS(chemicalSummary_conc_det) %>%
    mutate(guide_side = "Concentration [\U003BCg/L]") %>%
    left_join(cas_final, by="CAS"),
  graphData_tox_det = graph_chem_data_CAS(chemicalSummary_tox_det) %>%
    mutate(guide_side = "ToxCast [EAR]") %>%
    left_join(cas_final, by="CAS"),
  toxPlot_ear_conc_detects = fancy_combo(graphData_tox_det,
                                         graphData_conc_det,
                                         tox_list),
  chemicalSummary_conc_det_match = filter(chemicalSummary_conc_det, CAS %in% unique(graphData_tox_det$CAS)),
  graphData_conc_det_match = graph_chem_data_CAS(chemicalSummary_conc_det_match) %>%
    mutate(guide_side = "Concentration [\U003BCg/L]") %>%
    left_join(cas_final, by="CAS"),
  toxPlot_ear_conc_matches = fancy_combo(graphData_tox_det,
                                         graphData_conc_det_match,
                                         tox_list),
  graphData_conc_det_match_filter = graphData_conc_det_match,
  graphData_tox_det_filter = graphData_tox_det,
  toxPlot_ear_conc_matches_filter = fancy_combo(graphData_tox_det_filter,
                                                graphData_conc_det_match_filter,
                                                tox_list),
  AOP = readr::read_csv(file_in(readd(AOP_crosswalk))) %>%
                          select(endPoint=`Component Endpoint Name`, ID=`AOP #`) %>%
                          distinct(),
  aop_graph = plot_tox_endpoints_manuscript(chemicalSummary, AOP,
                                            category = "Chemical", 
                                            font_size = 7,title = " ",
                                            pallette = c("steelblue", "white")),
  site_info = prep_site_list(tox_list$chem_site),
  stack_plot = plot_tox_stacks_manuscript(chemicalSummary, 
                                          site_info, 
                                          font_size = 6,
                                          category = "Chemical Class")

)


# drake_config(data_analysis_plan)
# In R console:
# r_make("passive_data_setup.R")
# config <- drake_config(data_analysis_plan)
# vis_drake_graph(config, build_times = "none")
make(data_analysis_plan)

loadd(aop_graph)
pdf("plots/AOP.pdf", width = 4.5, height = 4.5)
ggarrange(
  aop_graph$count_plot, aop_graph$stackedPlot,
  aop_graph$chem_plot, aop_graph$aop_plot,nrow = 1,ncol = 4,
  widths =  c(2/10,6/10,1/10,1/10)
)
dev.off()

loadd(toxPlot_ear_conc)
pdf("plots/EAR_Conc_all.pdf", width = 4.5, height = 22, onefile=FALSE)
ggarrange(
  toxPlot_ear_conc$site_graph, 
  toxPlot_ear_conc$no_axis,
  widths =  c(3.25/9, 5.75/9),
  common.legend = TRUE, legend = "bottom"
)
dev.off()

loadd(toxPlot_ear_conc_detects)
pdf("plots/EAR_Conc_detects.pdf", width = 4.5, height = 11, onefile=FALSE)
ggarrange(
  toxPlot_ear_conc_detects$site_graph, 
  toxPlot_ear_conc_detects$no_axis,
  widths =  c(3.25/9, 5.75/9),
  common.legend = TRUE, legend = "bottom"
)
dev.off()

loadd(toxPlot_ear_conc_matches)
pdf("plots/EAR_Conc_detects_match.pdf", width = 4.75, height = 9, onefile=FALSE)
ggarrange(
  toxPlot_ear_conc_matches$site_graph, 
  toxPlot_ear_conc_matches$no_axis,
  widths =  c(3.25/9, 5.75/9),
  common.legend = TRUE, legend = "bottom"
)
dev.off()


###########################
loadd(graphData_tox_det)
loadd(graphData_conc_det_match)
loadd(tox_list)
source(file = "R/report/combo_plot2.R")
toxPlot_ear_conc_matches = fancy_combo(graphData_tox_det,
                                       graphData_conc_det_match,
                                       tox_list)

pdf("plots/EAR_Conc_detects_match.pdf", width = 4.5, height = 9, onefile=FALSE)
ggarrange(
  toxPlot_ear_conc_matches$site_graph, 
  toxPlot_ear_conc_matches$no_axis,
  widths =  c(3.25/9, 5.75/9),
  common.legend = TRUE, legend = "bottom"
)
dev.off()
