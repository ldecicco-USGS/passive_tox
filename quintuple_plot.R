
library(toxEval)
library(tidyverse)

path_to_data <- Sys.getenv("PASSIVE_PATH")


quad_fun <- function(path_to_data, tox_list){
  # source(file = here::here("read_chemicalSummary.R"))

  eco_list <- create_toxEval(file.path(path_to_data, "data/toxEval input file/","passive_benchmarks_all.xlsx"))
  
  tox_list$chem_info <- tox_list$chem_info %>% 
    rename(Chemical = chnm)
  
  chmn_df <-  tox_list$chem_info %>% 
    select(CAS, Chemical) %>% 
    distinct()
  
  eco_list$chem_info <- eco_list$chem_info %>% 
    rename(Chemical = chnm)
  
  summary_conc <- get_concentration_summary(tox_list)
  
  chmn_conc_df <- summary_conc %>% 
    select(CAS, chnm) %>% 
    distinct() %>%
    left_join(chmn_df, by = "CAS")
  
  orig_levels_conc <- data.frame(chnm = levels(chmn_conc_df$chnm)) %>% 
    left_join(chmn_conc_df %>% 
                mutate(chnm = as.character(chnm)), by = "chnm")
  
  cs_conc <- summary_conc %>% 
    select(-chnm) %>% 
    left_join(chmn_conc_df, by = "CAS") %>% 
    select(-chnm) %>%
    mutate(chnm = factor(Chemical, levels = orig_levels_conc$Chemical))
  
  gd_conc <- graph_chem_data(cs_conc)
  
  summary_eco <- get_chemical_summary(eco_list)
  
  chmn_eco_df <- summary_eco %>% 
    select(CAS, chnm) %>% 
    distinct() %>%
    left_join(chmn_df, by = "CAS")
  
  orig_levels_eco <- data.frame(chnm = levels(chmn_eco_df$chnm)) %>% 
    left_join(chmn_eco_df %>% 
                mutate(chnm = as.character(chnm)), by = "chnm")
  
  cs_eco <- summary_eco %>% 
    select(-chnm) %>% 
    left_join(chmn_eco_df, by = "CAS") %>% 
    select(-chnm) %>%
    mutate(chnm = factor(Chemical, levels = orig_levels_eco$Chemical))
  
  summary_eco_1 <- cs_eco %>% 
    filter(Bio_category == 1)
  
  summary_eco_2 <- cs_eco %>% 
    filter(Bio_category == 2)
  
  gd_eco_1 <- graph_chem_data(summary_eco_1)
  gd_eco_2 <- graph_chem_data(summary_eco_2)
  
  chmn_tox_df <-  chemicalSummary %>% 
    select(CAS, chnm) %>% 
    distinct() %>%
    left_join(chmn_df, by = "CAS")
  
  orig_levels_tox <- data.frame(chnm = levels(chmn_tox_df$chnm)) %>% 
    left_join(chmn_tox_df %>% 
                mutate(chnm = as.character(chnm)), by = "chnm")
  
  cs <- chemicalSummary %>% 
    select(-chnm) %>% 
    left_join(chmn_tox_df, by = "CAS") %>% 
    select(-chnm) %>%
    mutate(chnm = factor(Chemical, levels = orig_levels_tox$Chemical))
  
  gd_tox <- graph_chem_data(cs)
  
  combo1 <- side_by_side_data(gd_tox, gd_conc,
                              left_title = "ToxCast EAR",
                              right_title = "Concentration [\U003BCg/L]")
  
  
  combo2 <- side_by_side_data(gd_eco_1, gd_conc,
                             left_title = "ECOTOX Group 1 TQ",
                             right_title = "Concentration [\U003BCg/L]")
  
  combo3 <- side_by_side_data(gd_eco_2, gd_conc,
                              left_title = "ECOTOX Group 2 TQ",
                              right_title = "Concentration [\U003BCg/L]")
  
  combo_all <- combo1 %>%
    bind_rows(combo2 %>%
                filter(guide_side == "ECOTOX Group 1 TQ")
    ) %>%
    bind_rows(combo3 %>%
                filter(guide_side == "ECOTOX Group 2 TQ")
    )
  
  combo_all$Class <- factor(combo_all$Class,
                              levels = levels(combo1$Class))
  combo_all$chnm <- factor(combo_all$chnm,
                             levels = levels(combo1$chnm))
  combo_all$guide_side <- factor(combo_all$guide_side,
                                   levels = c("ToxCast EAR",
                                              "ECOTOX Group 1 TQ",
                                              "ECOTOX Group 2 TQ",
                                              "Concentration [\U003BCg/L]"))
  
  combo_all_no_NDs <- combo_all[combo_all$meanEAR > 0,]
  
  pretty_logs_with_check <- function(x){
    x <- x[!is.na(x)]
    pretty_range <- range(x[x > 0])
    pretty_logs <- 10^(-10:10)
    log_index <- which(pretty_logs < pretty_range[2] & pretty_logs > pretty_range[1])
    pretty_index_range <- diff(range(log_index))
    
    log_index <- c(log_index[1]-1, log_index, log_index[length(log_index)] + 1)
    pretty_logs_new <-  pretty_logs[log_index] 
  
    if(pretty_index_range > 6 & pretty_index_range < 9){
      pretty_logs_new <- pretty_logs_new[c(TRUE, FALSE)]
    } else if(pretty_index_range >= 9 & pretty_index_range < 12){
      pretty_logs_new <- pretty_logs_new[c(TRUE, FALSE, FALSE)]
    } else if(pretty_index_range >= 12){
      pretty_logs_new <- pretty_logs_new[c(TRUE, FALSE, FALSE, FALSE)]
    }
    
    return(pretty_logs_new)
  }
  
  class_colors <- function(chemicalSummary){
    
    # classes <- levels(chemicalSummary$Class)
    
    classes <- c("Insecticide","Flavor/Fragrance",
                 "Antimicrobial disinfectant", "Herbicide",
                 "Fire retardant", "Detergent metabolites",
                 "Pharmaceuticals", "Plasticizer",
                 "WW", "PAHs",
                 "Other", "OC Pesticides",
                 "Food Additive/Plasticizer", "Dye/Pigment",
                 "Solvent", "PBDEs", "Sterol", "Fuel", "PCBs")
    
    
    cbValues <- c("brown1", "gold","darkred","darkblue","yellow","grey10",
                  "darkolivegreen","darksalmon", "darkolivegreen1","cyan3","deeppink2","grey50",
                  "aquamarine","azure2","darkorange","darkorchid","cornflowerblue","cornsilk","white")
    
    names(cbValues) <- classes
    
    
    return(cbValues)
    
  }
  
  benchChems <- combo_all_no_NDs %>% 
    filter(as.character(guide_side) != "Concentration [\U003BCg/L]") %>% 
    select(chnm) %>% 
    distinct() %>% 
    mutate(chnm = as.character(chnm)) %>% 
    pull(chnm)
  
  combo_all_at_least2 <- combo_all_no_NDs 
  
    # filter(as.character(chnm) %in% benchChems) %>% 
    # mutate(chnm = droplevels(chnm))
    
  facet_labels <- combo_all_at_least2 %>% 
    select(guide_side) %>% 
    distinct() %>% 
    mutate(chnm = factor(tail(levels(combo_all_at_least2$chnm),1), levels = levels(combo_all_at_least2$chnm)),
           label = c("A", "D", "B", "C"),
           x = c(10^-7, 10^-6, 10^-9, 10^-9),
           thresh = c(0.001, NA, 0.1, 0.1))
  
  site_labels <- data.frame(
    guide_side = factor("ToxCast EAR", levels = levels(facet_labels$guide_side)),
    label = "# Sites",
    x = 10^-8,
    chnm = facet_labels$chnm[1]
  )
  
  counts <- combo_all_at_least2 %>% 
    group_by(chnm) %>% 
    summarize(n_sites = as.character(length(unique(site)))) %>% 
    ungroup() %>% 
    mutate(guide_side = factor("ToxCast EAR", levels = levels(facet_labels$guide_side)),
           x = 10^-8)

  test_plot4 <- ggplot() +
    geom_boxplot(data = combo_all_at_least2, 
                 aes(y = chnm, x = meanEAR, fill = Class),
                 lwd = 0.1, outlier.size = 0.1, na.rm = TRUE) +
    facet_grid(. ~ guide_side, scales = "free_x") +
    theme_minimal() +
    geom_vline(data = facet_labels,
               aes(xintercept = thresh), 
               colour="grey", linetype = "dashed",
               lwd = 0.5) +
    geom_text(data = facet_labels, 
              aes(y = chnm, label = label, x = x),
              hjust = "center", vjust = 1) +
    geom_text(data = site_labels, size = 2.5,
              aes(y = chnm, label = label, x = x),
              hjust = "center", vjust = -1) +
    geom_text(data = counts, size = 2,
              aes(y = chnm, label = n_sites, x = x),
              hjust = "center") +
    scale_fill_manual(values = class_colors(), drop=FALSE) +
    theme(axis.text.y = element_text(size = 6),
          axis.title = element_blank(),
          panel.grid.major = element_line(size = 0.1),
          legend.text =  element_text(size = 6),
          panel.border = element_rect(color = "black", fill = NA, size = 0.1),
          legend.margin = margin(t = 0, b = 0, 
                                 r = 0,l = 0)) +
    scale_x_log10(breaks = pretty_logs_with_check,
                  labels = toxEval:::fancyNumbers) +
    coord_cartesian(clip = "off") 
  
  return(test_plot4)

}

# ggsave(test_plot4, filename = "test_new_side.pdf", height = 11, width = 9)
# 
# 
# ggsave(test_plot4, filename = file.path(path_to_data,"Figures","quad_plot.pdf"),
#        height = 11, width = 9)


# Plot chemicals that aren't in either:
# 
# cs_conc_filtered <- cs_conc %>% 
#   filter(!(as.character(chnm) %in% benchChems),
#          EAR > 0)
# 
# cs_conc_filtered$Class <- droplevels(cs_conc_filtered$Class)
# 
# plot_chemical_boxplots(cs_conc_filtered, plot_ND = FALSE, 
#                        x_label =  "Concentration [\U003BCg/L]" )
# 
# library(formattable)
# library(DT)
# 
# # Table of chemical, class, and median:
# 
# fnc = function(var) {
#   var <- prettyNum(var)
#   var[var=="NA"] = ""
#   var
# }
# 
# library(sparkline)
# 
# gd_conc <- graph_chem_data(cs_conc) %>% 
#   filter(meanEAR > 0)
# 
# minEAR <- min(log10(gd_conc$meanEAR))
# maxEAR <- max(log10(gd_conc$meanEAR))
# 
# median_conc <- gd_conc %>% 
#   group_by(chnm, Class) %>% 
#   summarise(n_sites = length(unique(site)),
#             conc_median = median(meanEAR),
#             conc_box = spk_chr(log10(meanEAR),
#                                lineColor = 'black', 
#                                fillColor = '#ccc',
#                                chartRangeMin = minEAR,
#                                chartRangeMax = maxEAR,
#                                width = 200,
#                                height = 60,
#                                type = "box",
#                                #target = 8,
#                                highlightLineColor = 'orange', 
#                                highlightSpotColor = 'orange')) %>% 
#   ungroup() 
# 
# datatable(median_conc,
#           escape = F,
#           extensions = 'RowGroup',
#           options = list(rowGroup = list(dataSrc = 2),
#                          orderFixed = list(list(2,'asc')),
#                          pageLength = 100,
#                          fnDrawCallback = htmlwidgets::JS('function(){
#                                                           HTMLWidgets.staticRender();
#                                                           }'))) %>% 
#   spk_add_deps()
# 
# median_tox <- graph_chem_data(cs) %>% 
#   group_by(chnm, Class) %>% 
#   summarise(tox_median = median(meanEAR)) %>% 
#   ungroup() 
# 
# median_eco1 <- graph_chem_data(summary_eco_1) %>% 
#   group_by(chnm, Class) %>% 
#   summarise(eco1_median = median(meanEAR)) %>% 
#   ungroup() 
# 
# median_eco2 <- graph_chem_data(summary_eco_2) %>% 
#   group_by(chnm, Class) %>% 
#   summarise(eco2_median = median(meanEAR)) %>% 
#   ungroup() 
# 
# median_all <- median_conc %>% 
#   left_join(median_tox, by = c("chnm","Class")) %>% 
#   left_join(median_eco1, by = c("chnm","Class")) %>% 
#   left_join(median_eco2, by = c("chnm","Class")) %>% 
#   ungroup() %>% 
#   mutate(Class = factor(Class, levels = levels(cs_conc$Class)),
#          class_num = as.integer(Class),
#          conc_median = fnc(conc_median),
#          tox_median = fnc(tox_median),
#          eco1_median = fnc(eco1_median),
#          eco2_median = fnc(eco2_median))
# 
# 
# # mylog <- function(x){
# # 
# #   x <- as.numeric(x)
# # 
# #   x <- log10(x)
# #   
# #   x[!is.finite(x) & !is.na(x)] <- min(x[is.finite(x) & !is.na(x)])-2
# #   
# #   max_x = max(x, na.rm = TRUE)
# #   min_x = min(x, na.rm = TRUE)
# #   x <- (x - min_x)/(max_x - min_x)
# #   
# # }
# # 
# # 
# # df_f <- formattable(median_all,
# #             list(conc_median = color_bar("goldenrod", fun = mylog),
# #                  tox_median = color_bar("goldenrod", fun = mylog),
# #                  eco1_median = color_bar("goldenrod", fun = mylog),
# #                  eco2_median = color_bar("goldenrod", fun = mylog)))
# 
# styleColorBarLOG <- function(data, color, angle = 90){
# 
#   data <- log10(data)
#   rg = range(data, na.rm = TRUE, finite = TRUE)
#   r1 = rg[1]
#   r2 = rg[2]
#   r = r2 - r1
#   JS(sprintf("isNaN(Math.log10(parseFloat(value))) || Math.log10(value) <= %f ? '' : 'linear-gradient(%fdeg, transparent ' + (%f - Math.log10(value))/%f * 100 + '%%, %s ' + (%f - Math.log10(value))/%f * 100 + '%%)'", 
#              r1, angle, r2, r, color, r2, r))
# }
# 
# datatable(
#   median_all[,-5],
#   extensions = 'RowGroup',
#   options = list(rowGroup = list(dataSrc = 2),
#                  orderFixed = list(list(8,'asc')),
#                  pageLength = 100,
#                  autoWidth = TRUE,
#                  columnDefs = list(list(visible=FALSE, targets=c(2,8))))) %>% 
#   formatSignif(c(4:7), digits = 3) %>%  
#   formatStyle("conc_median",
#               background = styleColorBarLOG(as.numeric(median_all$conc_median), 'goldenrod'),
#               backgroundSize = '95% 80%',
#               backgroundRepeat = 'no-repeat',
#               backgroundPosition = 'center' ) %>% 
#   formatStyle("tox_median",
#               background = styleColorBarLOG(as.numeric(median_all$tox_median), 'goldenrod'),
#               backgroundSize = '95% 80%',
#               backgroundRepeat = 'no-repeat',
#               backgroundPosition = 'center' ) %>% 
#   formatStyle("eco1_median",
#               background = styleColorBarLOG(as.numeric(median_all$eco1_median), 'goldenrod'),
#               backgroundSize = '95% 80%',
#               backgroundRepeat = 'no-repeat',
#               backgroundPosition = 'center' ) %>%
#   formatStyle("eco2_median",
#               background = styleColorBarLOG(as.numeric(median_all$eco2_median), 'goldenrod'),
#               backgroundSize = '95% 80%',
#               backgroundRepeat = 'no-repeat',
#               backgroundPosition = 'center' )
# 
# library(sparkline)
