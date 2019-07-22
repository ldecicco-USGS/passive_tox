library(drake)
library(tidyverse)
library(googledrive)
library(readxl)
library(data.table)
library(toxEval)
library(openxlsx)

options(drake_make_menu = FALSE)

dir.create("data", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create(file.path("data","raw"), showWarnings = FALSE)
dir.create(file.path("data","clean"), showWarnings = FALSE)

# Set up Googledrive downloads:
source("R/setup/check_googledrive.R")
source("R/setup/file_config.R")
# Go from raw files to R objects:
source("R/analyze/data_reader.R")
source("R/analyze/get_sites_ready.R")
source("R/analyze/get_chem_info.R")
source("R/analyze/create_tox_file.R")

pkgconfig::set_config("drake::strings_in_dots" = "literals")

data_setup_plan <- drake_plan(
  pharm_2014_download = target(command = drive_download_gd(pharm_update_id,
                                              path = file_out("data/raw/pharm_update.xlsx"),
                                              time_stamp = last_modified_pharm),
                        trigger = trigger(change = last_modified_pharm)),
  chem_change_download = target(command = drive_download_gd(chem_name_id,
                                                           path = file_out("data/raw/chem_names.xlsx"),
                                                           time_stamp = last_modified_chem_name_id),
                               trigger = trigger(change = last_modified_pharm)),
  file_2014_download = target(command = drive_download_gd(general_2014_id,
                                                      path = file_out("data/raw/general_2014.xlsx"),
                                                      time_stamp = last_modified_general_2014),
                          trigger = trigger(change = last_modified_general_2014)),  
  file_2010_download = target(command = drive_download_gd(general_2010_id,
                                                          path = file_out("data/raw/general_2010.xlsx"),
                                                          time_stamp = last_modified_general_2010),
                              trigger = trigger(change = last_modified_general_2010)),
  WW_2014_download = target(command = drive_download_gd(ww_update_id,
                                                          path = file_out("data/raw/ww_update.xlsx"),
                                                          time_stamp = last_modified_ww_update),
                              trigger = trigger(change = last_modified_ww_update)),
  cas_download = target(command = drive_download_gd(cas_id,
                                                          path = file_out("data/raw/cas.xlsx"),
                                                          time_stamp = last_modified_cas_id),
                              trigger = trigger(change = last_modified_cas_id)), 
  cas_change_download = target(command = drive_download_gd(cas_change_id,
                                                    path = file_out("data/raw/cas_change.xlsx"),
                                                    time_stamp = last_modified_cas_change_id),
                        trigger = trigger(change = last_modified_cas_change_id)), 
  cas_change = readxl::read_xlsx(cas_change_download),
  cas_df = all_cas(cas_download),
  clean_cas_df = clean_cas(cas_df),
  clean_cas_fixed = fix_cas(clean_cas_df, cas_change),
  
  AOP_crosswalk = target(command = drive_download_gd(AOP_update_id,
                                                           path = file_out("data/raw/AOP_crosswalk.csv"),
                                                           time_stamp = last_modified_AOP),
                               trigger = trigger(change = last_modified_AOP)),
  site_download = target(command = drive_download_gd(site_id,
                                                    path = file_out("data/raw/sites_from_OWC.txt"),
                                                    time_stamp = last_modified_site_id),
                        trigger = trigger(change = last_modified_site_id)),  
  sites_OWC = data.table::fread(site_download,
                                data.table = FALSE, 
                                sep="\t", select = c("SiteID", "site_grouping", "Short Name"),
                                colClasses = c("SiteID"="character")),
  class_download = target(command = drive_download_gd(class_id,
                                                      path = file_out("data/raw/chem_classes.csv"),
                                                      time_stamp = last_modified_class_id),
                          trigger = trigger(change = last_modified_class_id)),  
  chem_info_old = read.csv(file_in("data/raw/chem_classes.csv"), stringsAsFactors = FALSE),
  exclude_download = target(command = drive_download_gd(exclude_id,
                                                        path = file_out("data/raw/exclude.csv"),
                                                        time_stamp = last_modified_exclude_id),
                            trigger = trigger(change = last_modified_exclude_id)),  
  exclude = get_exclude(file_in("data/raw/exclude.csv")),
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
  pharm_2014 = generic_file_opener(file_in("data/raw/pharm_update.xlsx"), cas_df,
                                    n_max = 41, 
                                    sheet = "est water concentrations",
                                    site_sheet = "PrioritySiteInfo",
                                    year = 2014,
                                    skip = 7, skip_site = 2),
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
  chem_info = get_chem_info(all_data, chem_info_old),
  all_data_fixed_cas = fix_cas(all_data, cas_change),
  chem_info_fixed_cas = fix_cas(chem_info, cas_change),
  saveRDS(object = chem_info_fixed_cas, file = file_out("data/clean/cas_df.rds")),
  sites = get_sites_ready(file_2014_download, file_2010_download, sites_OWC),
  tox_list_init = create_tox_object(all_data_fixed_cas, chem_info_fixed_cas, sites, exclude),
  saveOutput = openxlsx::write.xlsx(tox_list_init, file = file_out("data/clean/passive.xlsx"))
  
)

make(data_setup_plan)


config <- drake_config(data_setup_plan)
vis_drake_graph(config, build_times = "none")
# sankey_drake_graph(config)


