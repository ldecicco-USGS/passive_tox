# Mixtures
# Modifying this script from what Steve did for mixtures
# in general this script:

# 1. For each sample, sums EARs by endpoint
# 2. For each chemical, calculates percent contribution to each summed endpoint
# 3. Chooses a cutoff for percent contribution, and finds important mixtures.
# 4. Defines chemicals in important mixtures
sum_endpoints <- function(all_EARs, ear_cutoff = 0.001) {
  
  
  # note Steve limited each chemical to contributing to one endpoint by 
  # only using the max EAR val for each site-date-chem combination
  # I don't really know the rationale for this -- so am leaving out for now
  # I've incorporated a slightly different filter below: I used only the max summed endpointEAR
  # per site-date to represent sort of the "worst" mixture. Should touch base with Steve about this. 
  summed_EARs <- all_EARs %>%
    group_by(site, shortName, date, endPoint) %>%
    mutate(sum_ear_endpoint = sum(EAR)) %>%
    ungroup() %>%
    mutate(chem_mix_contribution = (EAR/sum_ear_endpoint)*100) %>% 
    filter(sum_ear_endpoint > ear_cutoff,
            chem_mix_contribution > 1)
  
  return(summed_EARs)
}

calc_contr_chems <- function(summed_EARs) {
  
  EAR_sum_endpoint <- summed_EARs %>%
    filter(EAR > 0) %>% 
    # arrange(CAS) %>% 
    group_by(site, shortName, date, endPoint, sum_ear_endpoint) %>%
    summarize(n_contr_chems = n(),
              contr_chems_lt = list(as.character(unique(chnm))),
              contr_cas_lt = list(as.character(unique(CAS))),
              contr_chems = paste(as.character(unique(chnm)), collapse = ","),
              contr_chems = paste(as.character(unique(CAS)), collapse = ","),
              max_individual_contr = max(EAR)) %>%
    ungroup() 
    # filter(n_contr_chems > 1)
  
  return(EAR_sum_endpoint)
}

calc_top_mixtures <- function(EAR_sum_endpoint, max_only = TRUE) {
  
  if (max_only) {
    top_mixtures <- EAR_sum_endpoint %>%
      group_by(site, shortName, date) %>%
      summarize(max_sum_ear_endpoint = max(sum_ear_endpoint),
                endPoint_top = endPoint[which.max(sum_ear_endpoint)],
                n_contr_chems = n_contr_chems[which.max(sum_ear_endpoint)],
                contr_chems = contr_chems[which.max(sum_ear_endpoint)][order(contr_chems[which.max(sum_ear_endpoint)])],
                max_individual_contr = max_individual_contr[which.max(sum_ear_endpoint)]) %>%
      mutate(prop_ind_contr = max_individual_contr/max_sum_ear_endpoint)
  } else {
    top_mixtures <- EAR_sum_endpoint %>%
      mutate(prop_ind_contr = max_individual_contr/sum_ear_endpoint,
             endPoint_top = endPoint,
             max_sum_ear_endpoint = sum_ear_endpoint)
  }
  
  return(top_mixtures)
}

top_mixes_fn <- function(contributing_chems, n_site_thresh) {
  
  chm_key <- contributing_chems %>% 
    select(contr_chems) %>% 
    distinct() %>% 
    left_join(unique(select(contributing_chems,
                            contr_chems,
                            contr_chems_lt)), 
              by="contr_chems")

  top_mixes <- contributing_chems %>% 
    group_by(contr_chems, endPoint) %>% 
    summarise(n_samples = n(),
              unique_sites = length(unique(site))) %>% 
    filter(unique_sites > {{n_site_thresh}}) %>% 
    arrange(desc(n_samples)) %>% 
    left_join(chm_key, by="contr_chems") %>% 
    mutate(contr_chems_st =
             paste(unique(unlist(contr_chems_lt)),
                   collapse = ",\n")) %>% 
    filter(!duplicated(contr_chems)) %>% 
    ungroup()

  return(top_mixes)
}

get_combos <- function(chnm, EAR, ear_cutoff){
  
  unique_chems <- unique(as.character(chnm))
  length_chems <- length(unique_chems)
  EAR <- setNames(EAR, unique_chems)
  n_combos <- sapply(seq_len(length_chems), function(x) factorial(length_chems)/(factorial(x)*factorial(length_chems-x)))
  n_sums <- cumsum(n_combos)
  df_tots <- data.frame(chems = rep(NA_character_, sum(n_combos)),
                        EARsum = rep(NA, sum(n_combos)), 
                        n_chems = rep(NA, sum(n_combos)),
                        stringsAsFactors = FALSE)
  
  for(n_chems in seq_len(length(unique_chems))){
    chems <- combn(unique_chems, n_chems, simplify = FALSE)
    y <- sapply(chems, function(x) sum(EAR[x]))
    chems_char <- sapply(chems, function(x) paste0(x, collapse = "|")  )
    if(n_chems == 1){
      df_tots$chems[1:n_sums[n_chems]] <- chems_char
      df_tots$EARsum[1:n_sums[n_chems]] <- y
      df_tots$n_chems[1:n_sums[n_chems]] <- n_chems
    } else {
      df_tots$chems[(n_sums[n_chems-1]+1):n_sums[n_chems]] <- chems_char
      df_tots$EARsum[(n_sums[n_chems-1]+1):n_sums[n_chems]] <- y
      df_tots$n_chems[(n_sums[n_chems-1]+1):n_sums[n_chems]] <- n_chems
    }
  }
  
  df_tots <- df_tots %>% 
    filter(EARsum > {{ear_cutoff}}) 
  
  chems <- df_tots$chems
  names(chems) <- df_tots$n_chems
  
  return(df_tots)
}

all_mixes_fn <- function(EAR_sum_endpoint, ear_cutoff) {
  
  df <- EAR_sum_endpoint %>% 
    ungroup() %>% 
    filter(EAR > 0) %>% 
    group_by(endPoint, site, shortName, date) %>% 
    summarize(chems = list(get_combos(chnm, EAR, {{ear_cutoff}}))) %>% 
    unnest(chems) %>% 
    ungroup()

  CASs <- EAR_sum_endpoint %>% 
    ungroup() %>% 
    filter(EAR > 0) %>% 
    group_by(endPoint, site, shortName, date) %>% 
    summarize(CASs = list(get_combos(CAS, EAR, {{ear_cutoff}}))) %>% 
    unnest(CASs) %>% 
    select(endPoint, CASs = chems, EARsum, n_chems, date, site, shortName) %>% 
    ungroup()

  df <- df %>% 
    left_join(CASs, by = c("endPoint", "EARsum", "n_chems", "date", "site", "shortName")) %>% 
    group_by(endPoint, chems, CASs, n_chems) %>% 
    summarize(n_samples = n(),
              n_sites = length(unique(site)),
              siteIDs = paste0(site, collapse = "|"),
              siteNames = paste0(shortName, collapse = "|")) %>% 
    ungroup() 

  return(df)
  
}

get_final_mixtures <- function(chemicalSummary,
                               EAR_thresh,
                               site_thresh_percent, tox_list){
  
  EARsum_endpoint <- sum_endpoints(chemicalSummary,
                                   ear_cutoff = EAR_thresh)
  
  cas_check <- tox_list$chem_data %>% 
    group_by(CAS) %>% 
    summarize(n_sites = length(unique(SiteID)))
  
  top_mixes <- all_mixes_fn(EARsum_endpoint, EAR_thresh) %>% 
    rowwise() %>% 
    mutate(total_sites = !!cas_check %>% 
             filter(CAS %in% unlist(strsplit(split = "\\|", x = CASs))) %>% 
             pull(n_sites) %>% 
             min()) %>% 
    ungroup() %>% 
    mutate(site_percent = 100*n_sites/total_sites) %>% 
    filter(site_percent >= !!site_thresh_percent)

  join_everything <- join_everything_fnx(chemicalSummary)
  
  mix_df <- clean_top_mixes(join_everything, 
                            chemicalSummary,
                            top_mixes, 
                            site_thresh) %>% 
    select(Assay = endPoint, 
           AOPs,
           Genes = genes,
           Pathways = pathways,
           Chemicals = chem_list,
           Class,
           CASs = CAS_list,
           `# Sites` = n_sites,
           siteIDs = siteIDs,
           siteNames = siteNames)
  
  mix_df$Chemicals <- rapply(mix_df$Chemicals, function(x)
    ifelse(x == "Di(2-ethylhexyl) phthalate", "DEHP", x), 
    how = "replace")
  
  mix_df$Chemicals <- rapply(mix_df$Chemicals, function(x)
    ifelse(x == "Tris(1,3-dichloro-2-propyl) phosphate", "TDCPP", x), 
    how = "replace")
  
  mix_df$Chemicals <- rapply(mix_df$Chemicals, function(x)
    ifelse(x == "Tris(2-butoxyethyl) phosphate", "TBEP", x), 
    how = "replace")
  
  mix_df$Chemicals <- rapply(mix_df$Chemicals, function(x)
    ifelse(x == "Tris(2-chloroisopropyl)phosphate", "TCPP", x), 
    how = "replace")
  
  mix_df$Chemicals <- rapply(mix_df$Chemicals, function(x)
    ifelse(x == "4-(1,1,3,3-Tetramethylbutyl)phenol", "p-tert-Octylphenol", x), 
    how = "replace")
  
  return(mix_df)
}

create_Excel_wb_mix <- function(df){
  
  max_chems <- sapply(df$Chemicals, function(x)length(x))
  max_AOPs <-  sapply(df$AOPs, function(x)length(x))
  max_Class <-  sapply(df$Class, function(x)length(x))
  
  max_chems <- pmax(max_chems, max_AOPs, max_Class)
  
  
  wb <- createWorkbook()
  addWorksheet(wb, "Mixtures")
  header_st <- createStyle(textDecoration = "Bold")
  writeData(wb = wb, sheet =  "Mixtures", 
            x = df, sep = " \n", headerStyle = header_st)
  setRowHeights(wb, "Mixtures", 
                rows = 1:nrow(df)+1, 
                heights = 15*max_chems)
  setColWidths(wb, "Mixtures", cols = 7, widths = 25)
  setColWidths(wb, "Mixtures", cols = 5, widths = 25)
  setColWidths(wb, "Mixtures", cols = 2, widths = 25)
  setColWidths(wb, "Mixtures", cols = c(1,3:4, 6), widths = "auto")
  addStyle(wb, sheet = "Mixtures", cols = 1:7, 
           rows = 1:nrow(df)+1, gridExpand = TRUE,
           style = createStyle(valign = 'top', wrapText = TRUE))
  return(wb)
}

class_key_fnx <- function(chemicalSummary){
  x <- chemicalSummary %>%
    select(chnm, Class, CAS) %>%
    distinct() %>%
    mutate(chnm = as.character(chnm),
           Class = as.character(Class),
           CAS = as.character(CAS))
  return(x)
}

clean_top_mixes <- function(join_everything, 
                            chemicalSummary,
                            top_mixes, 
                            n_sites){
  
  library(tidyverse)
  
  class_key <- class_key_fnx(chemicalSummary)
  
  df <- join_everything  %>% 
    left_join(top_mixes, by = "endPoint") %>% 
    group_by(endPoint, AOPs, genes, pathways) %>%
    filter(n_chems == max(n_chems)) %>%
    ungroup() %>% 
    mutate(chem_list = c(strsplit(chems, split = "\\|")),
           chem_list = sapply(chem_list, function(x) x[!(x %in% "")]),
           chem_list = sapply(chem_list, function(x) unique(x)),
           CAS_list = c(strsplit(CASs, split = "\\|")),
           CAS_list = sapply(CAS_list, function(x) x[!(x %in% "")]),
           CAS_list = sapply(CAS_list, function(x) unique(x)),
           AOPs = c(strsplit(AOPs, split = ",")),
           AOPs =  sapply(AOPs, function(x) x[!(x %in% "")]),
           AOPs = sapply(AOPs, function(x) unique(x)),
           max_n_chems = sapply(chem_list, function(x) length(x)))

  #Never do this! :)
  df$Class <- NA
  for(i in 1:nrow(df)){
    CAS_i <- df$CAS_list[i][[1]]
    df$Class[i] <- list(filter(class_key, CAS %in% CAS_i) %>% select(Class) %>% distinct() %>% pull(Class))
  }
  
  df$AOPs[sapply(df$AOPs, function(x) length(x) == 0)] <- ""
  
  return(df)
}

# calculate metrics by chemical
# calculate the number of times a chemical was in a mixture
# exclude 1-compound mixtures first

calc_chem_mix_metrics <- function(top_mixtures, summed_EARs, out_file) {
  
  top_mixes_2plus <- filter(top_mixtures, n_contr_chems >1)
  
  top_mix_chems <- summed_EARs %>%
    left_join(select(top_mixes_2plus, site, shortName, date, endPoint = endPoint_top, prop_ind_contr)) %>%
    filter(!is.na(prop_ind_contr)) %>%
    group_by(chnm) %>%
    summarize(times_in_mixes = n(),
              contribution_median = median(chem_mix_contribution),
              n_endpoints = length(unique(endPoint)),
              n_sites = length(unique(site))) %>%
    arrange(-times_in_mixes)
  
  write.csv(top_mix_chems, out_file, row.names = FALSE)
}

summarize_mixtures <- function(top_mixtures) {
  mix_summary <- top_mixtures %>%
    filter(n_contr_chems > 1) %>%
    group_by(contr_chems, n_contr_chems) %>%
    summarize(endPoint = paste(unique(endPoint_top), collapse = ', '), 
              mix_n_hits = n(),
              mix_n_hits_sites = length(unique(site)),
              mix_n_hits_samples = length(unique(paste0(site, date))),
              mix_n_hits_months = length(unique(lubridate::month(date))),
              mix_max_sum_EAR = max(max_sum_ear_endpoint),
              mix_median_sum_EAR = median(max_sum_ear_endpoint),
              mix_median_prop_individual_contr = median(prop_ind_contr))
  
}

summarize_by_n <- function(top_mixtures) {
  n_summary <-top_mixtures %>%
    filter(n_contr_chems >1) %>%
    group_by(n_contr_chems) %>%
    summarize(n_hits = n(),
              n_unique_mixes = length(unique(contr_chems)))
}

plot_mix_summary <- function(n_summary, mix_summary, top_mixtures, ear_sum, out_file) {
  p1 <- ggplot(n_summary, aes(x = n_contr_chems, y = n_hits)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label = n_unique_mixes), vjust = -0.25, color = 'red3') +
    scale_x_continuous(breaks = 2:8) +
    theme_bw() +
    labs(x = 'Number of chemicals in mixture', y = 'Hits (EARmix > 0.001)') +
    annotate('text', x = 6.6, y = 250, label = 'Number of \nunique mixtures', color = 'red3')
  
  top <- mix_summary %>% 
    filter(mix_n_hits_samples >= 10) %>%
    arrange(-mix_n_hits_samples)
  
  top_all <- filter(top_mixtures, contr_chems %in% top$contr_chems)  %>%
    mutate(top = TRUE)
  
  
  top_all_chems <- left_join(ear_sum, select(top_all, site, date, endPoint_top, contr_chems, top), by = c('site', 'date', 'endPoint' = 'endPoint_top')) %>%
    filter(!is.na(top))
  
  plot_summary <- top_all_chems %>%
    group_by(chnm, CAS, Class, contr_chems) %>%
    summarize(median_ear = median(EAR)) %>%
    left_join(top, by = 'contr_chems')
  
  plot_summary$contr_chems <- factor(plot_summary$contr_chems, levels = top$contr_chems)
  
  totals <- plot_summary %>%
    group_by(contr_chems, mix_n_hits_samples) %>%
    summarize(total = sum(median_ear)) %>%
    mutate(chnm = NA)
  
  plot_summary$chnm <- gsub('2,4-Dichlorophenoxyacetic acid', '2,4-D', plot_summary$chnm)
  p2 <- ggplot(plot_summary, aes(x = contr_chems, y = median_ear, fill = chnm)) +
    geom_bar(stat = 'identity', color = 'black') +
    geom_text(data = totals, aes(x = contr_chems, y = total, label = mix_n_hits_samples), vjust = -0.25, color = 'red3') +
    theme_bw() +
    geom_hline(yintercept = 0.001, linetype = 2) +
    theme(axis.text.x = element_blank(),
          legend.position = c(0.25, 0.7),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key.size = unit(0.3, "cm")) +
    labs(x = 'Unique mixtures', y = 'Median EARmix', fill = '') +
    annotate('text', x = 8, y = 0.016, label = 'Number of \nsample hits', color = 'red3')
  
  
  ggsave(out_file, cowplot::plot_grid(p1, p2), height = 3.6, width = 8)
  
  
}

# sites - calculate max EARmix across samples, as well as the 
# number of months where there is an EARmix > 0.001 
calc_site_mix_metrics <- function(top_mixtures, out_file) {
  site_mix <- top_mixtures %>%
    group_by(site, shortName) %>%
    summarize(n_mix_hits = n(),
              n_mix_hit_months = length(unique(lubridate::month(date))),
              max_EARmix = round(max(max_sum_ear_endpoint), 4)) %>%
    arrange(-max_EARmix)
  
  write.csv(site_mix, out_file, row.names = FALSE)
}


plot_trees <- function(form, sub_df, endpoint){
  
  set.seed(999)
  
  tree_2 <- rpart(formula = form,
                  data = sub_df,
                  control = rpart.control(minsplit = 10,
                                          minbucket = 10,
                                          cp = 0.001))
  
  tree_2.Prune <- prune(tree_2,
                        cp = tree_2$cptable[
                          which.min(tree_2$cptable[,"xerror"]), "CP"])
  
  plot(as.party(tree_2.Prune),
       tp_args = list(id=FALSE),
       main = paste(endpoint,"\n",
                    unique(sub_df$mix_st)))
  
  tree_terms <- unique(as.character(tree_2.Prune$frame$var)[as.character(tree_2.Prune$frame$var) != "<leaf>"])
  
  return(tree_terms)
}

get_lm_stuff <- function(form, sub_df, sumEAR, log){

  basic_lm <- lm(data = sub_df, formula = form)
  predictions <- predict(basic_lm, 
                         interval = 'confidence')
  
  predictions <- data.frame(predictions)
  
  sub_df_sub <- 
    sub_df %>% 
    bind_cols(predictions) %>% 
    select(sumEAR = {{sumEAR}}, #lwr, upr, #if we want a ribbon...bring this back 
           fit, shortName)
  
  if(log){
    sub_df_sub$fit <- 10^(sub_df_sub$fit)
    # sub_df_sub$lwr <- 10^(sub_df_sub$lwr)
    # sub_df_sub$upr <- 10^(sub_df_sub$upr)
  }
  
  x <- coef(basic_lm)
  x_df <- data.frame(x)
  x_df$coef <- row.names(x_df)
  
  dirty_eqn <- paste(format(x, digits = 2),
                     names(x), collapse = " +\n", sep = "*")
  dirty_eqn <- gsub("\\*\\(Intercept\\)", "", dirty_eqn)
  dirty_eqn <- gsub("+ -", "- ", dirty_eqn)
  
  if(log){
    dirty_eqn <- paste0("log10(",sumEAR, ") = ", dirty_eqn)
  } else {
    dirty_eqn <- paste0(sumEAR, " = ", dirty_eqn)
  }
  
  sub_df_sub <- sub_df_sub %>% 
    mutate(model = "lm")
  
  eqn_lm <- data.frame(x = 0,
                       y = 1.01*max(sub_df_sub$fit),
                       label = dirty_eqn,
                       model = "lm",
                       stringsAsFactors = FALSE)
  
  return(list(data_lm = sub_df_sub, label = eqn_lm, x_df = x_df))
  
}

get_surv_stuff <- function(form, sub_df, log){
  
  if(length(attr(terms(form),"term.labels")) == 0){
    return(list(data_lm = data.frame()), label = NA)
    message("No terms in survival")
  }
  
  basic_lm <- survival::survreg(form,
                                data = sub_df,
                                dist = "gaus")
  predictions <- predict(basic_lm)
  if(log){
    predictions <- data.frame(fit = 10^predictions)
  } else {
    predictions <- data.frame(fit = predictions)
  }
  
  sub_df_sub <- 
    sub_df %>% 
    bind_cols(predictions) %>% 
    select(lowEAR, highEAR, fit, shortName)
  
  x <- coef(basic_lm)
  x_df <- data.frame(x)
  x_df$coef <- row.names(x_df)
  
  dirty_eqn <- paste(format(x, digits = 2),
                     names(x), collapse = " +\n", sep = "*")
  dirty_eqn <- gsub("\\*\\(Intercept\\)", "", dirty_eqn)
  dirty_eqn <- gsub("+ -", "- ", dirty_eqn)
  
  if(log){
    dirty_eqn <- paste0("Surv(log10(lowEAR),log10(highEAR)) = ", dirty_eqn)
  } else {
    dirty_eqn <- paste0("Surv(lowEAR,highEAR) = ", dirty_eqn)
  } 
  
  sub_df_sub <- sub_df_sub %>% 
    mutate(model = "survival")
  
  eqn_surv <- data.frame(x = 0,
                       y = 1.01*max(sub_df_sub$fit),
                       label = dirty_eqn,
                       model = "survival",
                       stringsAsFactors = FALSE)
  
  return(list(data_surv = sub_df_sub, label = eqn_surv, x_df = x_df))

}

plot_lm <- function(form_lm, form_surv, sub_df, sumEAR = "sumEAR", 
                    log = FALSE){
  
  sub_df_lm_list <- get_lm_stuff(form = form_lm,
                            sub_df = sub_df, 
                            sumEAR = sumEAR, 
                            log = log)

  sub_df_surv_list <- get_surv_stuff(form = form_surv,
                                     sub_df = sub_df,
                                     log = log)
                      
  graph_data <- dplyr::bind_rows(sub_df_lm_list[["data_lm"]],
                                 sub_df_surv_list[["data_surv"]])

  graph_data$sumEAR[is.na(graph_data$sumEAR)] <- (graph_data$lowEAR[is.na(graph_data$sumEAR)] + graph_data$highEAR[is.na(graph_data$sumEAR)])/2
  
  label_info <- dplyr::bind_rows(sub_df_lm_list[["label"]],
                                 sub_df_surv_list[["label"]])
  
  x_df <- list(lm = sub_df_lm_list[["x_df"]],
               surv = sub_df_surv_list[["x_df"]])
  
  if(log){
    label_info$y = 10^((log10(min(graph_data$fit)) - log10(max(graph_data$fit)))/2)
  } else {
    label_info$y = 0.5*max(label_info$y)
  }
  
  
  scatter_urb_ag <- ggplot() +
    geom_point(data = graph_data,
               aes(x = sumEAR, y = fit)) +
    geom_text(data = filter(graph_data, sumEAR == max(sumEAR)),
              aes(x = sumEAR, y = fit, label = shortName),
              vjust = 1, hjust = 1) +
    geom_abline(slope = 1, 
                linetype = "solid", color = "red") +
    facet_wrap(. ~ model) +
    geom_text(data = label_info,
              aes(x = x, y = y, label = label),
              hjust = 0) +
    xlab("Observed") + ylab("Predicted") +
    theme_bw()

  if(log){
    scatter_urb_ag <- scatter_urb_ag +
      scale_y_log10() +
      scale_x_log10()
  }
    
  print(scatter_urb_ag)
  
  return(x_df)
}

get_formula <- function(sub_df, variables_to_use, 
                        sumEAR = "sumEAR", log = FALSE,
                        lasso = FALSE, survival = TRUE){
  
  if(lasso){
    predictors_df <- sub_df[,variables_to_use]
    response_df <- sub_df[,sumEAR]
    # from glmnet
    x <- as.matrix(predictors_df) # Removes class
    
    if(log){
      y <- log10(as.double(sub_df[[sumEAR]])) # Only class
      obs_term <- paste("log10(",sumEAR,")")
    } else {
      y <- as.double(sub_df[[sumEAR]]) # Only class
      obs_term <- sumEAR
    }
    
    # Fitting the model (Lasso: Alpha = 1)
    set.seed(999)
    cv.lasso <- cv.glmnet(x, y, #, family='multinomial',
                          alpha=1, 
                          standardize=TRUE, type.measure='mse')
    coefs_lasso <- as.matrix(coef(cv.lasso, s = "lambda.min"))
    coefs_to_save <- row.names(coefs_lasso)[which(coefs_lasso !=0)]
  
    bestlam = cv.lasso$lambda.min # Select lamda that minimizes training MSE
    
    lasso_coef = predict(cv.lasso,
                         type = "coefficients",
                         s = bestlam)
    coefs_to_save <- coefs_to_save[coefs_to_save != "(Intercept)"]
    
  } else {
    
    if(!survival){
      
      if(log){
        
        obs_term <- paste("log10(",sumEAR,")")
        
        form <- formula(paste(obs_term, " ~ ",
                              paste(variables_to_use,
                                    collapse = " + ")))
      } else {
        obs_term <- sumEAR
        
        form <- formula(paste(obs_term, " ~ ",
                              paste(variables_to_use,
                                    collapse = " + ")))
      }
      lm_high <- lm(form, data = sub_df)
      step_return <- MASS::stepAIC(lm_high, k = log(nrow(sub_df)), trace = FALSE)

    } else {
      
      if(log){
        obs_term <- 'survival::Surv(log10(lowEAR),
                             log10(highEAR), 
                             type="interval2")'        
      } else {
        obs_term <- 'survival::Surv(lowEAR, highEAR, 
                             type="interval2")'
      }

      
      form <- reformulate(termlabels = variables_to_use, 
                          response = obs_term)

      surv_high <- survival::survreg(form,
                                     data = sub_df,
                                     dist = "gaus")
      
      step_return <- MASS::stepAIC(surv_high, k = log(nrow(sub_df)), trace = FALSE)

    }
    
    coefs_to_save <- coefficients(step_return)
    coefs_to_save <- names(coefs_to_save)
    coefs_to_save <- coefs_to_save[coefs_to_save != "(Intercept)"]

  }
  
  if(length(coefs_to_save) == 0){

      new_form <- formula(paste(obs_term,"~ 1"))

  } else {
    new_form <- formula(paste0(obs_term, " ~ ",
                               paste(coefs_to_save, 
                                     collapse = " + ")))
  }
  return(new_form)
}

variable_summary <- function(chem, endpoint, x_df, x_df2){
  
  chems <- unlist(chem)
  lm_lin <- x_df[["lm"]]$coef[x_df[["lm"]]$coef != "(Intercept)"]
  lm_lin <- gsub("_", " ", lm_lin)

  lm_log <- x_df2[["lm"]]$coef[x_df2[["lm"]]$coef != "(Intercept)"]
  lm_log <- gsub("_", " ", lm_log)

  surv_lin <- x_df[["surv"]]$coef[x_df[["surv"]]$coef != "(Intercept)"]
  surv_lin <- gsub("_", " ", surv_lin)
  
  surv_log <- x_df2[["surv"]]$coef[x_df2[["surv"]]$coef != "(Intercept)"]
  surv_log <- gsub("_", " ", surv_log)
  
  if(length(lm_lin) == 0){
    lm_lin <- ""
  }
  if(length(lm_log) == 0){
    lm_log <- ""
  }
  if(length(surv_lin) == 0){
    surv_lin <- ""
  }
  if(length(surv_log) == 0){
    surv_log <- ""
  }
  
  df_temp <- data.frame(matrix("", ncol = 5, 
                               nrow = max(sapply(list(lm_lin,
                                                      lm_log,
                                                      surv_lin,
                                                      surv_log,
                                                      c(chems, endpoint)), length))
                               ),
                        stringsAsFactors = FALSE
                        )
  
  names(df_temp) <- c("Chem/Endpoint", "lm", "lm_log", "surv", "surv_log")
  
  
  df_temp$`Chem/Endpoint`[1:length(chems)] <- chems
  df_temp$`Chem/Endpoint`[length(chems) + 1] <- endpoint
  df_temp$lm[1:length(lm_lin)] <- lm_lin
  df_temp$lm_log[1:length(lm_log)] <- lm_log
  df_temp$surv[1:length(surv_lin)] <- surv_lin
  df_temp$surv_log[1:length(surv_log)] <- surv_log
  
  return(df_temp)
}


create_censor <- function(df){
  # TODO: actually us MDL EARS!!!
  offset <- 0.5*min(df$sumEAR[df$sumEAR > 0])
  
  df$lowEAR <- df$sumEAR + offset
  df$highEAR <- df$sumEAR + offset
  
  # Linear:
  df$highEAR[df$sumEAR == 0] <- 2*offset
  return(df)
}

