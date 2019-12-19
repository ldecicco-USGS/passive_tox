library(drake)
library(tidyverse)
library(toxEval)
library(ggpubr)

loadd(chemicalSummary)
loadd(site_info)
loadd(tox_list)
loadd(cas_final)

source(file = "R/report/stack_plots.R")
source(file = "R/report/combo_plot2.R")
color_map <- class_colors(tox_list)
font_size <- 8

loadd(graphData_tox_det)
ordered_class <- levels(graphData_tox_det$Class)

levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Food Additive/Plasticizer"] <- "Food Additive"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial disinfectant"] <- "Antimicrobial"

ordered_class[ordered_class == "Food Additive/Plasticizer"] <- "Food Additive"
ordered_class[ordered_class == "Antimicrobial disinfectant"] <- "Antimicrobial"

tox_list$chem_info$Class[tox_list$chem_info$Class == "Food Additive/Plasticizer"] <- "Food Additive"
tox_list$chem_info$Class[tox_list$chem_info$Class == "Antimicrobial disinfectant"] <- "Antimicrobial"

name_key <- chemicalSummary %>% 
  select(CAS, chnm) %>% 
  distinct() %>% 
  left_join(select(cas_final, CAS, chnm_new = chnm), by = "CAS") %>% 
  arrange(chnm)

name_key$chnm_new[name_key$chnm_new == "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"] <- "Tris(1,3-dichloro-2-propyl)phosphate"
name_key$chnm_new[name_key$chnm_new == "Tris(2-butoxyethyl)phosphate (TBEP)"] <- "Tris(2-butoxyethyl)phosphate"
name_key$chnm_new[name_key$chnm_new == "Tris(1-chloro-2-propyl)phosphate (TCPP)"] <- "Tris(1-chloro-2-propyl)phosphate"
name_key$chnm_new[name_key$chnm_new == "Tris(2-ethylhexyl)phosphate (TEHP)"] <- "Tris(2-ethylhexyl)phosphate"
                    
                    
levels(chemicalSummary$chnm) <- name_key$chnm_new[which(name_key$chnm == levels(chemicalSummary$chnm))]

classes <- unique(tox_list$chem_info$Class)

num_chem_to_keep <- 5

class_plots <- list()

site_info <- site_info %>% 
  arrange(`Short Name`)

levels(site_info$`Short Name`) <- paste0(levels(site_info$`Short Name`),
                                         " (",site_info$map_nm,")")

for(i in ordered_class){
  
  chem_i <- chemicalSummary %>%
    dplyr::filter(Class == i)
  
  if(nrow(chem_i) == 0){
    next
  }
  chem_i$chnm <- droplevels(chem_i$chnm)
  
  if(length(levels(chem_i$chnm)) > num_chem_to_keep){
    too_small <- levels(chem_i$chnm)[1:(length(levels(chem_i$chnm))-num_chem_to_keep)]
    just_right <- levels(chem_i$chnm)[(length(levels(chem_i$chnm))-num_chem_to_keep+1):length(levels(chem_i$chnm))]
    too_small_text <- paste0("Other [",length(too_small),"]")
    chem_i$chnm <- as.character(chem_i$chnm)
    chem_i$chnm[chem_i$chnm %in% too_small] <- too_small_text
    chem_i$chnm <- factor(chem_i$chnm, levels = c(rev(just_right), too_small_text))
    
  } else {
    chem_i$chnm <- factor(as.character(chem_i$chnm), levels = rev(levels(chem_i$chnm)))
  } 
  
  color_i <- colorspace::rainbow_hcl(length(unique(chem_i$chnm)), 
                                     start = -360, end = -55, 
                                     c = 100, l = 64)
  
  names(color_i) <- levels(chem_i$chnm)
  
  if(grepl("Other",names(color_i)[length(color_i)])){
    color_i[length(color_i)] <- "grey77"
  }
  
  class_plots[[i]] <- whole_stack(chem_i, site_info, tox_list,
                            color_i, font_size, 
                            category = "Chemical",
                            title = i)
}


library(cowplot)

general_lp <- c(0.05,0.87)

pdf("plots/MoreChemicalStacks_new.pdf", width = 10, height = 8, onefile = TRUE)

rel_widths <- c(3,3,rep(c(1,3),2), 1, 4)

plot_grid(
  class_plots[[ordered_class[1]]]$chem_count,
  class_plots[[ordered_class[1]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[ordered_class[2]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[2]]]$no_axis +
    theme(legend.position = general_lp ,
          strip.text.y =  element_blank()) +
    ylab("Sum of Maximum EAR Values"),
  class_plots[[ordered_class[3]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[3]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[ordered_class[4]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[4]]]$no_axis +
    theme(legend.position = c(0.015, 0.84),
          strip.text.y = element_text(size = 0.75*font_size)) +
    ylab(""),
  rel_widths = rel_widths,
  nrow = 1
)

plot_grid(
  class_plots[[ordered_class[5]]]$chem_count ,
  class_plots[[ordered_class[5]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[ordered_class[6]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[6]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab("Sum of Maximum EAR Values"),
  class_plots[[ordered_class[7]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[7]]]$no_axis +
    theme(legend.position = general_lp ,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[ordered_class[8]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[8]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y = element_text(size = 0.75*font_size))+
    ylab("") ,
  rel_widths = rel_widths,
  nrow = 1
)

plot_grid(
  class_plots[[ordered_class[9]]]$chem_count ,
  class_plots[[ordered_class[9]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[ordered_class[10]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[10]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab("Sum of Maximum EAR Values"),
  class_plots[[ordered_class[11]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[11]]]$no_axis +
    theme(legend.position = c(0.025, 0.92) ,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[ordered_class[12]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[12]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y = element_text(size = 0.75*font_size))+
    ylab("") ,
  rel_widths = rel_widths,
  nrow = 1
)

plot_grid(
  class_plots[[ordered_class[13]]]$chem_count ,
  class_plots[[ordered_class[13]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y = element_blank()) +
    ylab(""),
  class_plots[[ordered_class[14]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[14]]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y = element_blank()) +
    ylab("Sum of Maximum EAR Values"),
  class_plots[[ordered_class[15]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[15]]]$no_axis +
    theme(legend.position = c(0.05, 0.92) ,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[ordered_class[16]]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[ordered_class[16]]]$no_axis +
    theme(legend.position = c(0.05, 0.90),
          strip.text.y = element_text(size = 0.75*font_size))+
    ylab("") ,
  rel_widths = rel_widths,
  nrow = 1
)

dev.off()

