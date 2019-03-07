library(drake)
library(tidyverse)
library(googledrive)
library(readxl)
library(data.table)
library(toxEval)
library(openxlsx)

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
file_out_data <- file.path("data","clean","passive.xlsx")

data_setup_plan <- drake_plan(

  pharm_2014_download = target(command = drive_download_gd(pharm_update_id,
                                              path = file_out("data/raw/pharm_update.xlsx"),
                                              time_stamp = last_modified_pharm),
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
  cas_df = all_cas(cas_download),
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
  chem_info_old = read.csv(class_download, stringsAsFactors = FALSE),
  chem_info = get_chem_info(all_data, chem_info_old),
  exclude_download = target(command = drive_download_gd(exclude_id,
                                                        path = file_out("data/raw/exclude.csv"),
                                                        time_stamp = last_modified_exclude_id),
                            trigger = trigger(change = last_modified_exclude_id)),  
  exclude = get_exclude(exclude_download),
  OC_2014 = generic_file_opener(file_2014_download, cas_df,
                                 n_max = 45, 
                                 sheet = "OC-PCB-PBDE",
                                 site_sheet = "site info",
                                 year = 2014),
  PAHs_2014 = generic_file_opener(file_2014_download, cas_df,
                                   n_max = 33, 
                                   sheet = "PAHs",
                                   site_sheet = "site info",
                                   year = 2014),
  pharm_2014 = generic_file_opener(pharm_2014_download, cas_df,
                                    n_max = 41, 
                                    sheet = "est water concentrations",
                                    site_sheet = "PrioritySiteInfo",
                                    year = 2014,
                                    skip = 7, skip_site = 2),
  PAHs_2010 = generic_file_opener(file_2010_download, cas_df, 
                                   n_max = 33,
                                   sheet = "PAHs",
                                   site_sheet = "site info",
                                   year = 2010,
                                   skip_site = 2),
  OC_2010 = generic_file_opener(file_2010_download, cas_df,
                                 n_max = 40, 
                                 sheet = "OC-PCB-PBDE",
                                 site_sheet = "site info",
                                 year = 2010,
                                 skip_site = 2),
  WW_2010 = generic_file_opener(file_2010_download, cas_df, 
                                 n_max = 53, 
                                 sheet = "WW",
                                 site_sheet = "site info",
                                 year = 2010,
                                 skip_site = 2),
  WW_2014 = generic_file_opener(WW_2014_download, cas_df, 
                                 n_max = 46, 
                                 sheet = "est water concentrations",
                                 site_sheet = "PrioritySiteInfo",
                                 year = 2014,
                                 skip = 7,
                                 skip_site = 2),
  pharm_2010 = generic_file_opener(file_2010_download, cas_df,
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
  sites = get_sites_ready(file_2014_download, file_2010_download, sites_OWC),
  tox_list = create_tox_object(all_data, chem_info, sites, exclude),
  openxlsx::write.xlsx(tox_list, file = file_out(file_out_data), append=TRUE)
  
)

make(data_setup_plan)

config <- drake_config(data_setup_plan)
vis_drake_graph(config, build_times = "none")
# sankey_drake_graph(config)


