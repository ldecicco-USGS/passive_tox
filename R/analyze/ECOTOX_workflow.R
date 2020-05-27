# Ecotox workflow:
dir.create("R/analyze/plots", showWarnings = FALSE)
dir.create("R/analyze/out", showWarnings = FALSE)

#Combine ECOTOX endpoints, filter out for relevance, run stats
source("R/analyze/ecotox_combined_analyis.R")

source("R/analyze/ecotox_stats.R")


#graph TQs
#source("R/analyze/TQ_boxplots.R")
#source("R/analyze/graph_ECOTOX_by_effect_group.R")

#Analyze results: determine priorities
source("R/analyze/explore_threshold_exceedances_ECOTOX.R")
