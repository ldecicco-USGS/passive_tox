

#get_mixtures_bio_pathways_table <- function() {
library(here)
source(here("R/mixtures/mix_script.R"))
source(here("R/mixtures/prepare_mixture_data.R"))
chemicalSummary <- readRDS(file = file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","clean","chemical_summary.rds"))


n_sites <- 7 #10% of sites
EAR_thresh <- 0.001
EARsum_endpoint <- sum_endpoints(chemicalSummary,
                                 ear_cutoff = EAR_thresh)
top_mixes <- all_mixes_fn(EARsum_endpoint, EAR_thresh)
df <- clean_top_mixes(join_everything, 
                      top_mixes, 
                      n_sites)
#return(df)
#}



# chem_list <- character()
# for(i in 1:length(df$chem_list)) chem_list <- c(chem_list,df$chem_list[[i]])
# 
# unique(chem_list)
# 
# df$chem_list[[1]]
