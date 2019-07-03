loadd(graphData_conc_det_match)
loadd(graphData_tox_det)
loadd(chemicalSummary_conc)
loadd(tox_list)
loadd(cas_final)

# cas_final =  cas_df 
axis_num <- 6

source(file = "R/report/combo_plot2.R")

chemicalSummary_conc_no_match = chemicalSummary_conc %>%
  filter(!(CAS %in% unique(graphData_tox_det$CAS)),
         EAR > 0)

graphData_conc_no_match = graph_chem_data_CAS(chemicalSummary_conc_no_match) %>%
  mutate(guide_side = "Concentration [\U003BCg/L]") %>%
  left_join(cas_final, by="CAS")

full_classes <- c(levels(graphData_tox_det$Class),
                  levels(graphData_conc_no_match$Class)[!(levels(graphData_conc_no_match$Class) %in% levels(graphData_tox_det$Class))])

graphData_tox_det$Class <- factor(as.character(graphData_tox_det$Class), levels = full_classes)
graphData_conc_no_match$Class <- factor(as.character(graphData_conc_no_match$Class), levels = full_classes)

matches <- fancy_combo(graphData_tox_det, 
                       graphData_conc_det_match, 
                       tox_list, 
                       axis_size = axis_num)

# ggarrange(
#   matches$site_graph,matches$no_axis,
#   common.legend = TRUE, legend = "bottom"
# )

n_chems_matches <- length(unique(graphData_tox_det$chnm))

graphData_empty <- graphData_conc_no_match[FALSE,]

gd_no_match <- combine_gd(graphData_conc_no_match, graphData_empty)

n_chems_no_match <- length(unique(gd_no_match$chnm))

color_map <- class_colors(tox_list)
toxPlot_no_match <- combo_plot_matches_2(gd_no_match,
                                   axis_size = axis_num,
                                   color_map)
text_df_c <- label_info(gd_no_match, labels_to_use = "C")
toxPlot_no_match_w_lab <- add_label(toxPlot_no_match, text_df_c)
no_axis_no_match <- strip_graph(toxPlot_no_match_w_lab)
site_counts_df_no_match <- site_counts(tox_list$chem_data, no_axis_no_match$data)
site_graph_no_match <- site_count_plot(site_counts_df_no_match,
                                       axis_size = axis_num)

# ggarrange(
#   site_graph_no_match,no_axis_no_match,
#   common.legend = TRUE, legend = "bottom"
# )

# pdf("plots/triple_graph_half.pdf", width = 9, height = 5.5, onefile=FALSE)
# ggarrange(
#   
#   matches$site_graph,
#   matches$no_axis,
#   ggarrange(
#     site_graph_no_match, 
#     no_axis_no_match,
#     NULL,
#     nrow = 2, ncol = 2,
#     widths = c(1,1.4),
#     heights = c(n_chems_no_match,n_chems_matches-n_chems_no_match)
#   ),
#   widths =  c(1.5,4,4),nrow=1,ncol=3,
#   common.legend = TRUE, legend = "bottom"
# )
# dev.off()


pdf("plots/triple_graph_full_page.pdf", width = 9, height = 11, onefile=FALSE)
ggarrange(
  
  matches$site_graph,
  matches$no_axis,
  ggarrange(
    site_graph_no_match, 
    no_axis_no_match,
    NULL,
    nrow = 2, ncol = 2,
    widths = c(1.75,2),
    heights = c(n_chems_no_match,n_chems_matches-n_chems_no_match)
  ),
  widths =  c(2,4,4),nrow=1,ncol=3,
  common.legend = TRUE, legend = "bottom"
)
dev.off()

