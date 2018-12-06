library(toxEval)

tox_list_concentrations <- create_toxEval("cleanedData/passive_conc.xlsx")
chemicalSummary_conc <- get_chemical_summary(tox_list_concentrations)
