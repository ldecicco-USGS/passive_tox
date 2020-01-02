library(drake)
library(tidyverse)
library(toxEval)
library(ggpubr)

loadd(chemicalSummary)
loadd(site_info)
loadd(tox_list)
source(file = "R/report/stack_plots.R")
source(file = "R/report/combo_plot2.R")
source(file = "R/analyze/table_chem_class_Land_use_correlation.R")

lakes_ordered <- c("Lake Superior",
                   "Lake Michigan",
                   "Lake Huron",
                   "Lake Erie",
                   "Lake Ontario")

df <- Chem_Class_correlation_table()

df <- tibble(x = 0.95,
             y = 0.95,
             site_grouping = factor("Lake Superior", 
                                    levels = lakes_ordered),
             tb = list(df))

color_map <- class_colors(tox_list)
font_size <- 5

upperPlot <- plot_tox_stacks_manuscript2(chemical_summary = chemicalSummary, 
                                         chem_site = site_info, 
                                         title=NA,cbValues = color_map,
                                         font_size =  font_size, 
                                         category = "Chemical Class")

library(ggpmisc)

stack2 <- upperPlot +
  geom_table_npc(data = df,
                 aes(npcx = x, 
                     npcy = y,
                     label = tb), size = 2,
                 hjust = 1, vjust = 1)

ggsave(stack2, filename = "plots/new_stack_w_table.pdf", height = 8, width = 5)

full_plot <- whole_stack(chemicalSummary, site_info, title=NA,
                         tox_list, color_map, font_size, 
                         category = "Chemical Class")



pdf("plots/stack_full.pdf", width = 5.5, height =7, onefile=FALSE)
ggarrange(
  full_plot$chem_count,
  full_plot$no_axis +
    ylab("Sum of Maximum EAR Values"),
  widths = c(1.2,5),
  common.legend = TRUE, legend = "bottom"
)
dev.off()


levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Food Additive/Plasticizer"] <- "Food Additive"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial disinfectant"] <- "Antimicrobial"

tox_list$chem_info$Class[tox_list$chem_info$Class == "Food Additive/Plasticizer"] <- "Food Additive"
tox_list$chem_info$Class[tox_list$chem_info$Class == "Antimicrobial disinfectant"] <- "Antimicrobial"


classes <- unique(tox_list$chem_info$Class)
num_chem_to_keep <- 5
font_size <- 5
pdf("plots/classes_split_new.pdf", width = 6, height = 7, onefile=TRUE)
for(i in classes){

  chem_i <- dplyr::filter(chemicalSummary, Class == i)
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
                                     start = -360, end = -55, c = 100, l = 64)
  
  names(color_i) <- levels(chem_i$chnm)
  
  if(grepl("Other",names(color_i)[length(color_i)])){
    color_i[length(color_i)] <- "grey77"
  }
  
  class_plot <- whole_stack(chem_i, site_info, tox_list,
                            color_i, font_size, 
                            category = "Chemical",
                            title = i)
  
  i <- gsub("/","_",i)
  # pdf(paste0(i,".pdf"), width = 4.5, height = 5.5, onefile=FALSE)
  
  print(
  ggarrange(
    class_plot$chem_count,
    class_plot$no_axis +
      ylab("Sum of Maximum EAR Values"),
    widths = c(1.3,5),
    common.legend = TRUE, legend = "bottom"
  ))
  
}

dev.off()

class_plots <- list()

site_info <- site_info %>% arrange(`Short Name`)
levels(site_info$`Short Name`) <- paste0(levels(site_info$`Short Name`),
                                         " (",site_info$map_nm,")")

font_size <- 8

for(i in classes){
  
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

general_lp <- c(0.1,0.87)


pdf("plots/MoreChemicalStacks_new.pdf", width = 10, height = 8)

rel_widths <- c(3.75,3,rep(c(1,3),3), 1, 3.5)

plot_grid(
  class_plots[[1]]$chem_count,
  class_plots[[1]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[2]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[2]]$no_axis +
    theme(legend.position = c(0.1,0.92) ,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[3]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[3]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab("Sum of Maximum EAR Values"),
  class_plots[[4]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[4]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[5]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[5]]$no_axis +
    theme(legend.position = c(0.01,0.84),
          strip.text.y = element_text(size = 0.75*font_size)) +
    ylab(""),
  rel_widths = rel_widths,
  nrow = 1
)

plot_grid(
  class_plots[[6]]$chem_count,
  class_plots[[6]]$no_axis +
    theme(legend.position = c(0.05,0.87),
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[7]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[7]]$no_axis +
    theme(legend.position = general_lp ,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[8]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[8]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab("Sum of Maximum EAR Values"),
  class_plots[[9]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[9]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[10]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[10]]$no_axis +
    theme(legend.position = c(0.01,0.84),
          strip.text.y = element_text(size = 0.75*font_size)) +
    ylab(""),
  rel_widths = rel_widths,
  nrow = 1
)

plot_grid(
  class_plots[[11]]$chem_count,
  class_plots[[11]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[12]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[12]]$no_axis +
    theme(legend.position = general_lp ,
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[13]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[13]]$no_axis +
    theme(legend.position = general_lp,
          strip.text.y =  element_blank()) +
    ylab("Sum of Maximum EAR Values"),
  class_plots[[14]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[14]]$no_axis +
    theme(legend.position = c(0.01,0.92),
          strip.text.y =  element_blank()) +
    ylab(""),
  class_plots[[15]]$chem_count +
    theme(axis.text.y = element_blank()),
  class_plots[[15]]$no_axis +
    theme(legend.position = c(0.01,0.88),
          strip.text.y = element_text(size = 0.75*font_size)) +
    ylab(""),
  rel_widths = rel_widths,
  nrow = 1
)

dev.off()


# Mock up stack of concentrations that aren't in toxCast:

loadd(chemicalSummary)
loadd(chemicalSummary_conc)
loadd(site_info)
loadd(tox_list)

non_tox <- chemicalSummary_conc %>%
  filter(!(CAS %in% chemicalSummary$CAS))

source(file = "R/report/stack_plots.R")
source(file = "R/report/combo_plot2.R")
color_map <- class_colors(tox_list)
font_size <- 5

full_plot <- whole_stack(non_tox, site_info, title="Chemicals not in toxCast",
                         tox_list, color_map, font_size, 
                         category = "Chemical Class")
pdf("plots/conc_stack_nonTox.pdf", width = 5.5, height =7, onefile=FALSE)
ggarrange(
  full_plot$chem_count,
  full_plot$no_axis +
    ylab("Concentration [\U003BCg/L]"),
  widths = c(1.2,5),
  common.legend = TRUE, legend = "bottom"
)
dev.off()
