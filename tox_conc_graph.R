library(toxEval)
library(ggplot2)

source("get_data.R")
source("get_data_conc.R")

# Remove nondetects:
chemicalSummary <- chemicalSummary %>%
  filter(EAR > 0)

chemicalSummary_conc <- chemicalSummary_conc %>%
  filter(EAR > 0)

# Combo graph:
source("combo_graph_function.R")

graphData_tox <- graph_chem_data(chemicalSummary)
graphData_tox$guide_side <- "atop(ToxCast,Maximum~EAR[Chem]~per~Site)"

graphData_conc <- graph_chem_data(chemicalSummary_conc)
graphData_conc$guide_side <- "Concentration\n[μg/L]"

graphData_tox$guide_up <- "A"
graphData_conc$guide_up <- "A"

toxPlot_ear_conc <- combo_plot_matches(gd_1 = graphData_tox, 
                                 gd_2 = graphData_conc, 
                                 thres_1 = NA, 
                                 thres_2 = NA, 
                                 drop = TRUE)
