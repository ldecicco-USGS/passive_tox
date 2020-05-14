library(toxEval)

file_path <- file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_all.xlsx")

tox_list_bench <- create_toxEval(file_path)
summary_bench <- get_chemical_summary(tox_list_bench)

pdf(file.path("R","Analyze","plots","TQ_boxplots.pdf"),height = 11)
p <- plot_tox_boxplots(summary_bench, 
                  category = "Chemical", 
                  sum_logic = FALSE,
                  x_label = "Toxicity Quotient")
print(p)
dev.off()


path_to_file <- file.path(file.path(path_to_data, "data", "data_for_git_repo","clean", "passive.xlsx")) 
tox_list <- create_toxEval(path_to_file)
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC)

pdf(file.path("R","Analyze","plots","TQ_EAR_boxplots.pdf"),height = 11)
cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep)
summary_tox <- get_chemical_summary(tox_list, 
                                    ACC, 
                                    filtered_ep)
gd_tox <- graph_chem_data(summary_tox)
gd_bench <- graph_chem_data(summary_bench, 
                            sum_logic = FALSE)

combo <- side_by_side_data(gd_tox, gd_bench, 
                           left_title = "ToxCast", 
                           right_title = "ECOTOX")
p <- plot_chemical_boxplots(combo, guide_side,
                       x_label = "") +
  ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
print(p)
dev.off()
