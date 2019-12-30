library(here)
source(here("R/mixtures/mix_script.R"))
source(here("R/mixtures/prepare_mixture_data.R"))

drake::loadd(chemicalSummary)
n_sites <- 10
EAR_thresh <- 0.001

EARsum_endpoint <- sum_endpoints(chemicalSummary,
                                 ear_cutoff = EAR_thresh)

top_mixes <- all_mixes_fn(EARsum_endpoint, EAR_thresh)

df <- clean_top_mixes(join_everything, 
                      top_mixes, 
                      n_sites)
df_lu <- open_land_use()

big_enough <- function(x, thresh = 2){
  max(x, na.rm = TRUE) > thresh
}

df_lu_filtered <- df_lu %>% 
  select(-site) %>% 
  select_if(big_enough) %>% 
  bind_cols(select(df_lu, frac_2010, frac_2014))

M <- cor(df_lu_filtered)

lu_vars <- colnames(M)
lu_vars_rev <- rev(lu_vars)

exclude <- c()
exclude_rev <- c()

for(i in seq_along(lu_vars)){
  
  check_var <- lu_vars[i]
  
  if(check_var %in% exclude){
    next
  } else {
    check_cor <- names(which(M[,check_var] > 0.9 | 
                               M[,check_var] < -0.9))
    check_cor <- check_cor[check_cor != check_var]
    exclude <- c(exclude, check_cor)
  }
  
}

df_lu_filtered <- df_lu_filtered[,-which(names(df_lu_filtered) %in% exclude)]

M2 <- cor(df_lu_filtered)
df_lu_filtered <- bind_cols(df_lu[,c("site","Urban","Crops")],
                            df_lu_filtered) %>% 
  select(-Basin_Area_mi2, -Population_Density)

all_chemicals_in_mixtures <- unique(unlist(df$chem_list))

big_enough <- function(x, thresh = 10){
  max(x, na.rm = TRUE) > thresh
}

cs_mix <- chemicalSummary %>% 
  filter(chnm %in% all_chemicals_in_mixtures,
         EAR > 0) %>% 
  select(site) %>% 
  distinct() %>% 
  left_join(df_lu_filtered, by = "site") %>% 
  select(-site) %>% 
  select_if(big_enough)

auto_categories <- names(cs_mix)

for(i in seq_len(nrow(df))){
  
  chems <- unlist(df$chem_list[i])
  mixture <- paste(chems, collapse = ",")
  endpoint <- df$endPoint[i] 
  
  cat("\n")
  cat("\n##", mixture,"\n")
  
  mixture <- paste(chems, collapse = ",\n")
  
  sub_df <- chemicalSummary %>% 
    filter(endPoint == {{endpoint}},
           EAR > 0) %>% 
    group_by(site, shortName, date) %>% 
    summarize(sumEAR = sum(EAR)) %>% 
    group_by(site, shortName) %>% 
    summarize(maxEAR = max(sumEAR)) %>% 
    left_join(df_lu_filtered, by = "site") %>% 
    mutate(mix_st = mixture) %>% 
    ungroup()
  
}