loadd(chemicalSummary)
loadd(site_info)
loadd(tox_list)
source(file = "R/report/stack_plots.R")
source(file = "R/report/combo_plot2.R")
color_map <- class_colors(tox_list)
axis_num <- 5
plot_back <- plot_tox_stacks_manuscript(chemicalSummary,
                                        cbValues = color_map, 
                                        chem_site = site_info,
                                        font_size = axis_num,
                                        category = "Chemical Class")
no_axis_plot_back <- strip_graph(plot_back)
chem_df <- chem_counts(chemicalSummary, site_info)
chem_count_graph <- chem_count_plot(chem_df,
                                    axis_size = axis_num)

pdf("stack.pdf", width = 4.5, height = 5.5, onefile=FALSE)
ggarrange(
  chem_count_graph,
  no_axis_plot_back,widths = c(1.3,5),
  common.legend = TRUE, legend = "bottom"
)
dev.off()
