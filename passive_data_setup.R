library(drake)
library(tidyverse)
library(googledrive)
library(readxl)
library(data.table)
library(toxEval)
library(openxlsx)

path_to_data <- Sys.getenv("PASSIVE_PATH")

options(drake_make_menu = FALSE)

dir.create("data", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create(file.path("data","raw"), showWarnings = FALSE)
dir.create(file.path("data","clean"), showWarnings = FALSE)

# Go from raw files to R objects:
source(file = "R/analyze/data_reader.R")
source(file = "R/analyze/get_sites_ready.R")
source(file = "R/analyze/get_chem_info.R")
source(file = "R/analyze/create_tox_file.R")
source(file = "R/analyze/open_land_use.R")

pkgconfig::set_config("drake::strings_in_dots" = "literals")

data_setup_plan <- drake_plan(

  cas_change_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/cas_change.xlsx")),
                              file_out("data/raw/cas_change.xlsx")),
  cas_change = readxl::read_xlsx(file_in("data/raw/cas_change.xlsx")),

  chem_info_old_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/chem_classes.csv")),
                                 file_out("data/raw/chem_classes.csv")),
  chem_info_old = read.csv(file_in("data/raw/chem_classes.csv"), stringsAsFactors = FALSE),
 
  
  exclude_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/exclude.csv")),
                           file_out("data/raw/exclude.csv")),
  exclude = get_exclude(file_in("data/raw/exclude.csv")),
  
  AOP_crosswalk_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/AOP_crosswalk.csv")),
                                 file_out("data/raw/AOP_crosswalk.csv")),
  OC_2014_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/general_2014.xlsx")),
            file_out("data/raw/general_2014.xlsx")),
  OC_2014 = generic_file_opener(file_in("data/raw/general_2014.xlsx"), cas_df,
                                 n_max = 45, 
                                 sheet = "OC-PCB-PBDE",
                                 site_sheet = "site info",
                                 year = 2014),

  PAHs_2014 = generic_file_opener(file_in("data/raw/general_2014.xlsx"), cas_df,
                                   n_max = 33, 
                                   sheet = "PAHs",
                                   site_sheet = "site info",
                                   year = 2014),
  
  pharm_2014_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/pharm_update.xlsx")),
                              file_out("data/raw/pharm_update.xlsx")),
  pharm_2014 = generic_file_opener(file_in("data/raw/pharm_update.xlsx"), cas_df,
                                    n_max = 41, 
                                    sheet = "est water concentrations",
                                    site_sheet = "PrioritySiteInfo",
                                    year = 2014,
                                    skip = 7, skip_site = 2),
  
  PAHs_2010_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/general_2010.xlsx")),
                             file_out("data/raw/general_2010.xlsx")),
  PAHs_2010 = generic_file_opener(file_in("data/raw/general_2010.xlsx"), cas_df, 
                                   n_max = 33,
                                   sheet = "PAHs",
                                   site_sheet = "site info",
                                   year = 2010,
                                   skip_site = 2),
  OC_2010 = generic_file_opener(file_in("data/raw/general_2010.xlsx"), cas_df,
                                 n_max = 40, 
                                 sheet = "OC-PCB-PBDE",
                                 site_sheet = "site info",
                                 year = 2010,
                                 skip_site = 2),
  WW_2010 = generic_file_opener(file_in("data/raw/general_2010.xlsx"), cas_df, 
                                 n_max = 53, 
                                 sheet = "WW",
                                 site_sheet = "site info",
                                 year = 2010,
                                 skip_site = 2),
  
  WW_2014_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/ww_update.xlsx")),
                           file_out("data/raw/ww_update.xlsx")),
  WW_2014 = generic_file_opener(file_in("data/raw/ww_update.xlsx"), cas_df, 
                                 n_max = 46, 
                                 sheet = "est water concentrations",
                                 site_sheet = "PrioritySiteInfo",
                                 year = 2014,
                                 skip = 7,
                                 skip_site = 2),
  pharm_2010 = generic_file_opener(file_in("data/raw/general_2010.xlsx"), cas_df,
                                    n_max = 44, 
                                    sheet = "pharms",
                                    site_sheet = "site info",
                                    year = 2010,
                                    skip_site = 2),
  all_data = bind_rows(pharm_2010,
                        WW_2010,
                        OC_2010,
                        PAHs_2010,
                        pharm_2014,
                        WW_2014,
                        OC_2014,
                        PAHs_2014) ,
  cas_df_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/cas.xlsx")),
            file_out("data/raw/cas.xlsx")),
  cas_df = all_cas(file_in("data/raw/cas.xlsx")),
  clean_cas_df = clean_cas(cas_df),
  clean_cas_fixed = fix_cas(clean_cas_df, cas_change),
  all_data_chnm = clean_names(all_data),
  all_data_fixed_cas = fix_cas(all_data_chnm, cas_change),
  chem_info = get_chem_info(all_data_fixed_cas, chem_info_old),
  chem_info_fixed_cas = fix_cas(chem_info, cas_change),
  out_cas = saveRDS(object = chem_info_fixed_cas, 
          file = file_out("data/clean/cas_df.rds")),
  out_cas_sync = saveRDS(object = chem_info_fixed_cas, 
          file = file_out(!!file.path(path_to_data,"data/data_for_git_repo/clean/cas_df.rds"))),  
  sites_orig_2014 = readxl::read_excel(file_in("data/raw/general_2014.xlsx"),
                     sheet = "site info",
                     skip = 3),
  sites_OWC_copy = file.copy(file_in(!!file.path(path_to_data,"data/data_for_git_repo/raw/sites_from_OWC.txt")),
                             file_out(file_out("data/raw/sites_from_OWC.txt"))),
  sites_OWC = data.table::fread(file_in("data/raw/sites_from_OWC.txt"),
                                data.table = FALSE, 
                                sep="\t", select = c("SiteID", "site_grouping", "Short Name"),
                                colClasses = c("SiteID"="character")),
  sites_orig_2010 = readxl::read_excel(file_in("data/raw/general_2010.xlsx"),
                                       sheet = "site info",
                                       skip = 2),
  df_lu = open_land_use(),
  sites = get_sites_ready(sites_orig_2014, sites_orig_2010, sites_OWC, df_lu),
  tox_list_init = create_tox_object(all_data_fixed_cas, chem_info_fixed_cas, sites, exclude),

  saveOutput = openxlsx::write.xlsx(tox_list_init, file = file_out("data/clean/passive.xlsx")),
  saveOutput2 = openxlsx::write.xlsx(tox_list_init, 
                                     file = file_out(!!file.path(path_to_data,"data/data_for_git_repo/clean/passive.xlsx")))
  
)

config <- drake_config(data_setup_plan)
vis_drake_graph(config, build_times = "none")

make(data_setup_plan)



# sankey_drake_graph(config)


