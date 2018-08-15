library(toxEval)
tox_list <- create_toxEval("cleanedData/passive.xlsx")

ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info = end_point_info)
filtered_ep <- filter_groups(cleaned_ep)

lakes_ordered <- c("Lake Superior",
                   "Lake Michigan",
                   "Lake Huron",
                   "Lake Erie",
                   "Lake Ontario")

tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping,
                                           levels=lakes_ordered)

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

plot_tox_boxplots(chemicalSummary, category = "Chemical")
