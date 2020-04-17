update_table_2 <- function(path_to_data){
  wb <- loadWorkbook(file.path(path_to_data, "Supplemental", "Table S2 - Chemical list+conc ranges.xlsx"))
  
  df <- readWorkbook(wb, startRow = 5, check.names = TRUE)
  
  df_joined <- df %>% 
    select(Class, Analyte, CAS) %>% 
    filter(Analyte != "") %>% 
    full_join(tox_list$chem_info, by = "CAS") %>% 
    select(CAS, Class.x, Class.y)
  
  changes = data.frame(
    chems = c("Dibenzo[a,h]anthracene", "Bupropion", "Citalopram",
              "Duloxetine", "Methadone", "	Propranolol",
              "Sertraline", "tramadol", "	Tris(1-chloro-2-propyl)phosphate (TCPP)"),
    Dave = c("53-07-3", "34911-55-2", "59729-33-8",
             "116539-59-4", "76-99-3", "525-66-6",
             "79617-96-2", "27203-92-5", "26248-87-3"),
    Laura = c("53-70-3", "34841-39-9", "219861-08-2",
              "136434-34-9", "1095-90-5", "318-98-9",
              "79559-97-0", "36282-47-0", "13674-84-5"),
    stringsAsFactors = FALSE)
  
  # write.csv(changes, "change_cas.csv", row.names = FALSE)
  
  df <- df %>% 
    filter(!is.na(Analyte))
  
  new_CAS <- df %>% 
    select(CAS) %>% 
    left_join(changes, by=c("CAS"="Dave")) %>% 
    mutate(CAS = ifelse(!is.na(Laura), Laura, CAS)) %>% 
    filter(!is.na(CAS)) 
  
  orig_cas <- df$CAS
  
  df$CAS <- new_CAS$CAS
  
  # writeData(wb, sheet = "Table S2",
  #           startRow = 6, startCol = 3, 
  #           x = df$CAS)
  
  # changedCAS <- createStyle(fgFill = "steelblue2")
  # 
  # addStyle(wb, sheet = "Table S2",
  #          style = changedCAS,
  #          cols = 3, gridExpand = FALSE, 
  #          rows = 5 + which(!orig_cas == new_CAS))
  
  
  ALL_TOX_DATA <- readRDS(file.path(Sys.getenv("PASSIVE_PATH"),
                                    "data","data_for_git_repo","raw",
                                    "all_tox.rds"))
  
  chem_info <- tox_list$chem_info
  
  ALL_TOX_DATA_in_study <-  ALL_TOX_DATA %>% 
    select(CAS = casn, endPoint=aenm, modl_acc, flags, hitc) %>% 
    filter(CAS %in% chem_info$CAS) %>% 
    group_by(CAS) %>% 
    summarize(`Total ToxCast assays` = length(unique(endPoint)),
              `Assays with hits` = length(unique(endPoint[hitc == 1])))
  
  assays_left <- chemicalSummary %>% 
    select(CAS, endPoint) %>% 
    distinct() %>% 
    group_by(CAS) %>% 
    summarize(`Assays in study` = length(unique(endPoint)))
  
  chem_info <- chem_info %>% 
    left_join(ALL_TOX_DATA_in_study, by="CAS", ) %>% 
    left_join(assays_left, by="CAS") %>% 
    mutate(`Total ToxCast assays` = ifelse(is.na(`Total ToxCast assays`), 0, `Total ToxCast assays`),
           `Assays with hits` = ifelse(is.na(`Assays with hits`) & 
                                         `Total ToxCast assays` != 0, 0, 
                                       `Assays with hits`),
           `Assays in study` = ifelse(is.na(`Assays in study`) & 
                                        `Assays with hits` != 0, 0, 
                                      `Assays in study`))
  
  format_2 <- function(x, nd_text = "ND"){
    x_txt <- ifelse(!is.finite(x) | x == 0,
                    nd_text,
                    formatC(signif(x, digits = 2), digits=2, format = "fg")
    )
    
    x_txt <- gsub(" ", "", x_txt)
    return(x_txt)
  }
  
  df_tox <- df %>% 
    select(CAS, Analyte) %>% 
    left_join(chem_info, by="CAS") %>% 
    select(-Analyte, -chnm, -Class) %>% 
    mutate(`2010_MDL` = format_2(1000*`2010_MDL`, nd_text = "--"),
           `2010_MQL` = format_2(1000*`2010_MQL`, nd_text = "--"),
           `2014_MDL` = format_2(1000*`2014_MDL`, nd_text = "--"),
           `2014_MQL` = format_2(1000*`2014_MQL`, nd_text = "--"))
  
  
  # writeData(wb, sheet = "Table S2",
  #           startRow = 5, startCol = 24, 
  #           x = df_tox)
  
  chem_stats <- tox_list$chem_data %>% 
    mutate(Value = 1000 * Value) %>%
    group_by(CAS) %>% 
    mutate(tots_mean = format_2(mean(Value[Value > 0])),
           tots_median = format_2(median(Value[Value > 0]))) %>% 
    group_by(CAS, `Sample Date`, tots_mean, tots_median) %>% 
    summarise(min = min(Value[Value > 0]),
              max = max(Value),
              mean = mean(Value[Value > 0]),
              median = median(Value[Value > 0]),
              n_dets = sum(is.na(comment)),
              samples = n()) %>% 
    ungroup() %>% 
    mutate(max = format_2(max),
           min = format_2(min),
           mean = format_2(mean),
           median = format_2(median),
           n_dets = format_2(n_dets, nd_text = "0")) %>% 
    pivot_wider(id_cols = c("CAS", "tots_mean", "tots_median"), 
                names_from = `Sample Date`, values_fill = list(mean = "--",
                                                               min = "--",
                                                               max = "--",
                                                               median = "--",
                                                               samples = "0",
                                                               n_dets = "--"),
                values_from = c("mean", "min", "max", "median", "n_dets", "samples")) %>% 
    arrange(match(CAS, df$CAS)) %>% 
    left_join(df_tox, by = "CAS")  %>% 
    rowwise() %>% 
    mutate(min_tots = format_2(min(c(as.numeric(min_2010), as.numeric(min_2014)), na.rm = TRUE)),
           max_tots = format_2(max(c(as.numeric(max_2010), as.numeric(max_2014)), na.rm = TRUE))) %>% 
    select(CAS,
           `Total ToxCast assays`, `Assays with hits`, `Assays in study`,
           min_2010, max_2010, median_2010, mean_2010, n_dets_2010, samples_2010, `2010_MDL`, `2010_MQL`,
           min_2014, max_2014, median_2014, mean_2014, n_dets_2014, samples_2014, `2014_MDL`, `2014_MQL`,
           min_tots, max_tots, tots_median, tots_mean, sites_det, sites_tested)
  
  writeData(wb, sheet = "Table S2",
            startRow = 6, startCol =6, colNames = FALSE,
            x = select(chem_stats, -CAS))
  
  saveWorkbook(wb, file = file.path(path_to_data, "Supplemental", "Table2_Combo.xlsx"), overwrite = TRUE)
  
}