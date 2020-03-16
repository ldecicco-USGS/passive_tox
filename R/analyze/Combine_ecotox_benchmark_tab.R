library(openxlsx)
library(readxl)
library(tidyverse)

path_to_data <- Sys.getenv("PASSIVE_PATH")

benchmarks <- read_xlsx(file.path(path_to_data,  "data", "toxEval input file", "passive_benchmarks_chems_in_toxcast.xlsx"),sheet = 5)
benchmarks_non <- read_xlsx(file.path(path_to_data,  "data", "toxEval input file", "passive_benchmarks_non_toxcast.xlsx"),sheet = 5)
benchmark_tab <- full_join(benchmarks,benchmarks_non)


wb <- loadWorkbook(file.path(path_to_data, "data", "data_for_git_repo","clean", "passive.xlsx"))
addWorksheet(wb,sheetName = "Benchmarks")
writeData(wb,sheet = "Benchmarks",x=benchmark_tab)
saveWorkbook(wb,file=file.path(path_to_data, "data", "toxEval input file", "passive_benchmarks_ECOTOX.xlsx"),overwrite = TRUE)


