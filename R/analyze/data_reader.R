generic_file_opener <- function(file_name, cas_df, n_max, sheet, site_sheet,
                                year, skip = 6, skip_site = 3){
  
  data_wide <- read_excel(file_name,
                          sheet = sheet,
                          skip = skip, n_max = n_max)
  
  site_stuff <- read_excel(file_name,
                           sheet = site_sheet,
                           skip = skip_site)
  
  if("CERC Site #" %in% names(site_stuff)){
    site_stuff <- rename(site_stuff, SiteID=`CERC Site #`)
  } else {
    site_stuff <- site_stuff %>%
      rename(SiteID = `CERC ID`,
             STAID = `USGS Station ID`)
  }
  
  site_stuff$SiteID <- gsub("site ","",site_stuff$SiteID, ignore.case = TRUE)
  
  units <- names(data_wide)[-1:-2]
  if(isTRUE(sum(grepl("pg/L", units)) == length(units))){
    convert <- 1000000
  } else if (isTRUE(sum(grepl("ng/L", units)) == length(units))){
    convert <- 1000
  } else if (isTRUE(sum(grepl("ug/L", units)) == length(units))){
    convert <- 1
  } else {
    stop("Check units!")
  }
  
  if(sheet != "est water concentrations"){
    names_wide <- read_excel(file_name,
                             sheet = sheet,
                             skip = 3, n_max = 1) 
    names(data_wide)[4:length(names(names_wide))] <- names(names_wide)[4:length(names(names_wide))]
    names(data_wide)[1] <- "chnm"
    names(data_wide)[2] <- "MDL"
    names(data_wide)[3] <- "MQL"
    data_long <- data_wide %>%
      tidyr::gather(SiteID, Value, -chnm, -MDL, -MQL) 
  } else {
    names_wide <- read_excel(file_name,
                             sheet = sheet,
                             skip = 4, n_max = 1)
    sites_start <- which(names(names_wide) == "site 1")
    names(data_wide)[sites_start:length(names(names_wide))] <- names(names_wide)[sites_start:length(names(names_wide))]
    names(data_wide)[1] <- "chnm"
    
    if(sites_start == 4){
      names(data_wide)[2] <- "DL"
      names(data_wide)[3] <- "RL"
      data_long <- data_wide %>%
        tidyr::gather(SiteID, Value, -chnm, -DL,-RL)
      sheet <- "WW"
    } else {
      names(data_wide)[2] <- "RL"
      data_long <- data_wide %>%
        tidyr::gather(SiteID, Value, -chnm, -RL)
      sheet <- "pharms"
    }
  }
  data_long <- na.omit(data_long)
  data_long$comment <- ""
  data_long$comment[grep("<",data_long$Value)] <- "<"
  data_long$comment[grep("DNQ",data_long$Value)] <- "DNQ"
  data_long$Value <- gsub("DNQ","",data_long$Value)
  data_long$Value <- gsub("<","",data_long$Value)
  data_long$Value <- gsub("a","",data_long$Value)
  data_long$Value <- gsub("b","",data_long$Value)
  data_long$Value <- gsub("c","",data_long$Value)
  data_long$Value <- gsub(" ","",data_long$Value)
  data_long <- data_long[data_long$Value != "lostinfield",]
  data_long <- data_long[data_long$Value != "-----",]
  data_long <- data_long[data_long$Value != "'-----",]
  data_long$comment[which(data_long$Value == "ND")] <- "<"
  data_long$Value[which(data_long$Value == "ND")] <- data_long$MDL[which(data_long$Value == "ND")]
  data_long <- data_long[data_long$Value != "NA",]
  
  data_long$Value <- as.numeric(data_long$Value) 
  data_long$Value <- data_long$Value/convert
  data_long$generic_class <- sheet
  data_long$`Sample Date` <- year
  data_long$SiteID <- gsub("site ","",data_long$SiteID, ignore.case = TRUE)
  
  data_long <- filter(data_long, 
                      !(is.na(Value) & comment == ""))
  
  data_long <- data_long %>%
    mutate(chnm = tolower(chnm)) %>%
    left_join(cas_df, by="chnm") %>%
    mutate(chnm = tools::toTitleCase(chnm)) %>%
    left_join(select(site_stuff, SiteID, STAID), by="SiteID") 
  
  data_long$STAID[nchar(data_long$STAID) == 8 & substring(data_long$STAID, first = 1, last = 1) != "0"] <- dataRetrieval::zeroPad(data_long$STAID[nchar(data_long$STAID) == 8 & substring(data_long$STAID, first = 1, last = 1) != "0"], 9)
  data_long$STAID <- dataRetrieval::zeroPad(data_long$STAID, 8)
  
  data_long <- data_long %>%
    select(-SiteID) %>%
    rename(SiteID=STAID) %>%
    filter(!is.na(chnm),
           CAS != "---" | is.na(CAS),
           CAS != "-" | is.na(CAS))
  
  ignore_totals <- c("Total PCBs","Total Pcbs in Mg/l","Total Oc Pesticides")
  data_long <- data_long[!(data_long$chnm %in% ignore_totals),]
  
  data_long$CAS[data_long$chnm == "Buproprion"] <- "34841-39-9"
  data_long$CAS[data_long$chnm == "Nadolol"] <- "42200-33-9"
  data_long$chnm[data_long$chnm == "Tris(1,3-Dichloro-2-Propyl)Phosphate (t"] <- "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"
  data_long$CAS[data_long$chnm == "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"] <- "13674-87-8"
  
  data_long <- data_long[!(data_long$chnm %in% c("Tcpp_isomer","Tcpp Isomer")),]
  
  data_long$CAS[data_long$chnm == "Omeprazole + Esomprazole"] <- "73590-58-6"
  
  if(any(is.na(data_long$CAS))){
    message("Some CAS didn't match up")
  }
  
  # Get rid of censored data:
  data_long$Value[data_long$comment != ""] <- 0
  
  data_long$SiteID[data_long$SiteID == "04085790"] <- "04085721"
  
  return(data_long)
}

all_cas <- function(file_cas="raw/cas.xlsx"){
  
  cas_data <- data.frame()
  
  for(i in c("OC-PCB-PBDE","PAHs","CERC WW","LC8240","LC8069","More")){
    tab_i <- read_excel(file_cas, sheet = i, skip = 3)
    cas_data <- bind_rows(cas_data, tab_i)
  }
  
  cas_data <- cas_data[!is.na(cas_data$Analyte),]
  
  cas_data_cleaned <- select(cas_data, chnm=Analyte, CAS=`CAS Number`) %>%
    mutate(chnm = tolower(chnm)) %>%
    filter(!(CAS %in% c("-","---"))) %>%
    filter(!duplicated(chnm)) %>%
    arrange(chnm)
    
  last_row <- nrow(cas_data_cleaned)
  cas_data_cleaned[last_row+1,] <- c("dl-menthol","89-78-1")
  
  
  return(cas_data_cleaned)
}

clean_cas <- function(cas_df){
  
  cas_final =  cas_df %>%
    filter(!duplicated(CAS)) %>%
    mutate(chnm = tools::toTitleCase(chnm))
  
  cas_final$chnm[cas_final$chnm == "Deet"] <- "DEET"
  cas_final$chnm[cas_final$chnm == "Tcep"] <- "TCEP"
  cas_final$chnm[cas_final$chnm == "Tbep"] <- "TBEP"
  cas_final$chnm[cas_final$chnm == "Tdcpp"] <- "TDCPP"
  cas_final$chnm[cas_final$chnm == "Total Pcbs"] <- "Total PCBS"
  cas_final$chnm[cas_final$chnm == "O,p'-Ddd"] <- "o,p'-DDD"
  cas_final$chnm[cas_final$chnm == "P,p'-Ddd"] <- "p,p'-DDD"
  cas_final$chnm[cas_final$chnm == "Pentachloroanisole (Pca)"] <- "PCA"
  cas_final$chnm[cas_final$chnm == "Tributyl Phosphate (Tbp)"] <- "TBP"
  cas_final$chnm[cas_final$chnm == "Hydrochlorothiazide (Hctz)"] <- "HCTZ"
  cas_final$chnm[cas_final$chnm == "Tris(2−Chloroethyl)Phosphate (Tcep)"] <- "TCEP"
  cas_final$chnm[cas_final$chnm == "O,p'-Ddt"] <- "o,p'-DDT"
  cas_final$chnm[cas_final$chnm == "O,p'-Ddt"] <- "o,p'-DDT"
  cas_final$chnm[cas_final$chnm == "P,p'-Dde"] <- "p,p'-DDE"
  cas_final$chnm[cas_final$chnm == "P,p'-Ddt"] <- "p,p'-DDT"
  cas_final$chnm[cas_final$chnm == "O,p'-Dde"] <- "o,p'-DDE"
  cas_final$chnm[cas_final$chnm == "Tris(1-Chloro-2-Propyl)Phosphate (Tcpp)"] <- "TCPP"
  cas_final$chnm[cas_final$chnm == "Hexachlorobenzene (Hcb)"] <- "HCB"
  cas_final$chnm[cas_final$CAS == "77-93-0"] <- "Triethyl Citrate "
  cas_final$chnm[cas_final$CAS == "30306-93-5"] <- "Ethyl Citrate"
  cas_final$chnm[grep("Pbde-", cas_final$chnm)] <- gsub(pattern = "Pbde-",
                                                        replacement = "PBDE-",
                                                        cas_final$chnm[grep("Pbde-", cas_final$chnm)])
  cas_final <- rbind(cas_final, data.frame(CAS="34841-39-9",
                                           chnm="Bupropion",
                                           stringsAsFactors = FALSE))
  cas_final$chnm[cas_final$CAS == "34911-55-2"] <- "Bupropion hydrochloride"
  
  cas_final$chnm[grep(pattern = "Delta-Benzenehexachloride",cas_final$chnm)] <- "delta-Bhc"
  cas_final$chnm[grep(pattern = "Beta-Benzenehexachloride",cas_final$chnm)] <- "beta-Bhc"
  cas_final$chnm[grep(pattern = "Alpha-Benzenehexachloride", cas_final$chnm)] <- "alpha-Bhc"
  
  return(cas_final)
}