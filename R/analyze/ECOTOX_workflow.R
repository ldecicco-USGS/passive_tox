# Ecotox workflow:
dir.create("R/analyze/plots", showWarnings = FALSE)
dir.create("R/analyze/out", showWarnings = FALSE)

# Chemicals included in ToxCast
source("R/analyze/ecotox_toxcast.R")
#source("R/analyze/graph_toxcast_ECOTOX.R")


# Chemicals not included in ToxCast
source("R/analyze/ecotox_non_toxcast.R")
#source("R/analyze/graph_non_toxcast_ECOTOX.R")

#Combine ToxCast and non ToxCast results
source("R/analyze/Combine_ecotox.R")
source("R/analyze/Combine_ecotox_benchmark_tab.R")

#graph TQs
source("R/analyze/TQ_boxplots.R")
