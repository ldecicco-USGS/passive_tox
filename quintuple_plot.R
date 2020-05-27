library(toxEval)
library(tidyverse)

path_to_data <- Sys.getenv("PASSIVE_PATH")
source(file = "read_chemicalSummary.R")

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
                            left_title = "ToxCast",
                            right_title = "Conc")


combo2 <- side_by_side_data(gd_eco_1, gd_conc,
                           left_title = "Eco_1",
                           right_title = "Conc")

combo3 <- side_by_side_data(gd_eco_2, gd_conc,
                            left_title = "Eco_2",
                            right_title = "Conc")

combo_all <- combo1 %>%
  bind_rows(combo2 %>%
              filter(guide_side == "Eco_1")
  ) %>%
  bind_rows(combo3 %>%
              filter(guide_side == "Eco_2")
  )

combo_all$Class <- factor(combo_all$Class,
                            levels = levels(combo1$Class))
combo_all$chnm <- factor(combo_all$chnm,
                           levels = levels(combo1$chnm))
combo_all$guide_side <- factor(combo_all$guide_side,
                                 levels = c("ToxCast",
                                            "Eco_1",
                                            "Eco_2",
                                            "Conc"))



# test_plot <- plot_chemical_boxplots(combo_all, guide_side,
#                        x_label = "", plot_ND = FALSE) +
#   ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
# 
# ggsave(test_plot, filename = "test.pdf", height = 11, width = 9)



combo_all_no_NDs <- combo_all[combo_all$meanEAR > 0,]

pretty_logs_new <- toxEval:::prettyLogs(combo_all_no_NDs$meanEAR[!is.na(combo_all_no_NDs$meanEAR)])

pretty_logs_new <- pretty_logs_new[c(TRUE, FALSE, FALSE)]

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

# test_plot2 <- ggplot() +
#   geom_boxplot(data = combo_all_no_NDs, 
#                aes(y = chnm, x = meanEAR, fill = Class),
#                lwd = 0.1, outlier.size = 0.1, na.rm = TRUE) +
#   facet_grid(. ~ guide_side, scales = "free_x") +
#   theme_minimal() +
#   scale_fill_manual(values = class_colors(), drop=FALSE) +
#   theme(axis.text.y = element_text(size = 5),
#         axis.title = element_blank(),
#         panel.grid.major = element_line(size = 0.1),
#         legend.text =  element_text(size = 8)) +
#   scale_x_log10(breaks = pretty_logs_new, labels = toxEval:::fancyNumbers)
# 
# ggsave(test_plot2, filename = "test2.pdf", height = 9, width = 11)


benchChems <- combo_all_no_NDs %>% 
  filter(guide_side != "Conc") %>% 
  select(chnm) %>% 
  distinct() %>% 
  mutate(chnm = as.character(chnm)) %>% 
  pull(chnm)

combo_all_at_least2 <- combo_all_no_NDs %>% 
  filter(as.character(chnm) %in% benchChems) %>% 
  mutate(chnm = droplevels(chnm))
  

test_plot3 <- ggplot() +
  geom_boxplot(data = combo_all_at_least2, 
               aes(y = chnm, x = meanEAR, fill = Class),
               lwd = 0.1, outlier.size = 0.1, na.rm = TRUE) +
  facet_grid(. ~ guide_side, scales = "free_x") +
  theme_minimal() +
  scale_fill_manual(values = class_colors(), drop=FALSE) +
  theme(axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        legend.text =  element_text(size = 8)) +
  scale_x_log10(breaks = pretty_logs_new, labels = toxEval:::fancyNumbers)

ggsave(test_plot3, filename = "test5.pdf", height = 11, width = 9)


combo_eco <- bind_rows(mutate(gd_eco_1,
                              guide_side = "1"),
                       mutate(gd_eco_2,
                              guide_side = "2"))

test_plotECO <- ggplot() +
  geom_boxplot(data = combo_eco, 
               aes(y = chnm, x = meanEAR, fill = Class),
               lwd = 0.1, outlier.size = 0.1, na.rm = TRUE) +
  facet_grid(. ~ guide_side, scales = "free_x") +
  theme_minimal() +
  scale_fill_manual(values = class_colors(), drop=FALSE) +
  theme(axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        legend.text =  element_text(size = 8)) +
  scale_x_log10(breaks = pretty_logs_new, labels = toxEval:::fancyNumbers)
