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
    filter(sum_ear_endpoint > ear_cutoff) %>%  
    filter(chem_mix_contribution > 1)
  
  return(summed_EARs)
}

calc_contr_chems <- function(summed_EARs) {
  
  EAR_sum_endpoint <- summed_EARs %>%
    arrange(CAS) %>% 
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

plot_lm <- function(form, sub_df, sumEAR = "sumEAR", log=FALSE){
  
  if(length(attr(terms(form),"term.labels")) == 0){
    
    scatter_new <- ggplot() +
      geom_point(data = sub_df,
                 aes(x = {{sumEAR}}, y = Urban)) +
      theme_bw()
    print(scatter_new)
    break
  }
  
  basic_lm <- lm(data = sub_df, formula = form)
  
  predictions <- predict(basic_lm, 
                         interval = 'confidence')
  
  predictions <- data.frame(predictions)
  
  sub_df_sub <- 
    sub_df %>% 
    bind_cols(predictions) %>% 
    select(sumEAR = {{sumEAR}}, lwr, upr, fit, shortName)
  
  x <- coef(basic_lm)
  x_df <- data.frame(x)
  x_df$coef <- row.names(x_df)
  
  dirty_eqn <- paste(format(x, digits = 2),
                     names(x), collapse = " +", sep = "*")
  dirty_eqn <- gsub("\\*\\(Intercept\\)", "", dirty_eqn)
  dirty_eqn <- gsub("+ -", "- ", dirty_eqn)
  
  if(log){
    dirty_eqn <- paste0("log10(",sumEAR, ") = ", dirty_eqn)
  } else {
    dirty_eqn <- paste0(sumEAR, " = ", dirty_eqn)
  }

  if(log){
    
    sub_df_sub$fit <- 10^(sub_df_sub$fit)
    sub_df_sub$lwr <- 10^(sub_df_sub$lwr)
    sub_df_sub$upr <- 10^(sub_df_sub$upr)
    
  } 
  
  scatter_urb_ag <- ggplot() +
    geom_point(data = sub_df_sub,
               aes(x = sumEAR, y = fit)) +
    geom_text(data = filter(sub_df_sub, sumEAR == max(sumEAR)),
              aes(x = sumEAR, y = fit, label = shortName),
              vjust = 1, hjust = 1) +
    geom_ribbon(data = sub_df_sub,
                aes(x = sumEAR,
                    ymin = lwr,
                    ymax = upr),
                alpha=0.3) +
    geom_abline(slope = 1, 
                linetype = "solid", color = "red") +
    geom_text(aes(label = dirty_eqn,
                  x = 0.0,
                  y = max(sub_df_sub$upr)*1.01),
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

get_formula <- function(sub_df, variables_to_use, sumEAR = "sumEAR", log=FALSE){
  
  predictors_df <- sub_df[,variables_to_use]
  response_df <- sub_df[,sumEAR]
  
  # from glmnet
  x <- as.matrix(predictors_df) # Removes class
  
  if(log){
    y <- log10(as.double(sub_df[[sumEAR]])) # Only class
  } else {
    y <- as.double(sub_df[[sumEAR]]) # Only class
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
  
  if(length(coefs_to_save) == 0){
    
    return(NULL)
    
  } else {
  
    if(log){
      new_form <- formula(paste0("log10(",sumEAR, ") ~ ", paste(coefs_to_save, collapse = " + ")))
    } else {
      new_form <- formula(paste0(sumEAR, " ~ ", paste(coefs_to_save, collapse = " + ")))
    }
    
    return(new_form)
  }
}

variable_summary <- function(chem, endpoint, x_df, x_df2){
  df_temp <- data.frame(chems_endpoint = paste(
    paste(unlist(chem),
          collapse = ",\n"),endpoint, sep=",\n"),
    tree_return = paste(tree_return,
                        collapse = ",\n"),
    linear_return = paste(x_df$coef[x_df$coef != "(Intercept)"],
                          collapse = ",\n"),
    log_return = paste(x_df2$coef[x_df2$coef != "(Intercept)"],
                       collapse = ",\n"),
    stringsAsFactors = FALSE)
  
  df_temp$log_return <- gsub("_"," ", df_temp$log_return)
  df_temp$linear_return <- gsub("_"," ", df_temp$linear_return)
  df_temp$tree_return <- gsub("_"," ", df_temp$tree_return)
  
  return(df_temp)
}
