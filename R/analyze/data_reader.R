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
  
  #Add replicate flag?
  
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
  data_long <- dplyr::filter(data_long, !is.na(Value))
  data_long <- dplyr::filter(data_long, !is.na(chnm))
  data_long$comment <- ""
  data_long$comment[grep("<",data_long$Value)] <- "<"
  data_long$comment[grep("DNQ",data_long$Value)] <- "DNQ"
  data_long$Value <- gsub("DNQ","",data_long$Value)
  data_long$Value <- gsub("<","",data_long$Value)
  data_long$Value <- gsub("a","",data_long$Value)
  data_long$Value <- gsub("b","",data_long$Value)
  data_long$Value <- gsub("c","",data_long$Value)
  data_long$Value <- gsub(" ","",data_long$Value)
  data_long$Value <- gsub("E","",data_long$Value)
  data_long <- data_long[data_long$Value != "lostinfield",]
  data_long <- data_long[data_long$Value != "-----",]
  data_long <- data_long[data_long$Value != "'-----",]
  data_long <- data_long[data_long$Value != "nosmple",]
  data_long$comment[which(data_long$Value == "ND")] <- "<"
  data_long$Value[which(data_long$Value == "ND")] <- data_long$MDL[which(data_long$Value == "ND")]
  data_long <- filter(data_long, Value != "NA")
  data_long$Value[data_long$Value == "7.0000000000000007-2"] <- "0.7"
  
  data_long$Value <- as.numeric(data_long$Value) 
  data_long$Value <- data_long$Value/convert
  data_long$generic_class <- sheet
  data_long$`Sample Date` <- year
  data_long$SiteID <- gsub("site ","",data_long$SiteID, ignore.case = TRUE)
  
  # Premature taking out censored values?
  data_long <- filter(data_long,
                      !(is.na(Value) & comment == ""))
  
  # Dealing with 3 relicates!
  # data_wide$Date[(data_wide$SiteID %in% c("04086000","04119400","04208000")) &&
  #                  data_wide$Date == 2014]
  
  if(!("Station shortname" %in% names(site_stuff))){
    site_stuff$`Station shortname` <- "temp"
    site_stuff$`Station shortname`[grepl("Replicate",site_stuff$`Site ID`)] <- "Replicate"
  }
  
  data_long <- data_long %>%
    mutate(chnm = tolower(chnm)) %>%
    left_join(cas_df, by="chnm") %>%
    mutate(chnm = tools::toTitleCase(chnm)) %>%
    left_join(select(site_stuff, SiteID, STAID, `Station shortname`), by="SiteID") 
  
  #Keep replicates:
  # data_long$`Sample Date`[grepl(pattern = "Replicate",
  #                               x = data_long$`Station shortname`)] <-  data_long$`Sample Date`[grepl(pattern = "Replicate",
  #                                                                                                      x = data_long$`Station shortname`)] + 0.5 #This essentially makes a 2nd "sample" for the data
  #Lose replicates:
  data_long <- data_long %>% 
    filter(!grepl(pattern = "Replicate",
                 x = `Station shortname`))
  
  data_long$STAID[nchar(data_long$STAID) == 8 & substring(data_long$STAID, first = 1, last = 1) != "0"] <- dataRetrieval::zeroPad(data_long$STAID[nchar(data_long$STAID) == 8 & substring(data_long$STAID, first = 1, last = 1) != "0"], 9)
  data_long$STAID <- dataRetrieval::zeroPad(data_long$STAID, 8)
  
  data_long <- data_long %>%
    select(-SiteID, -`Station shortname`) %>%
    rename(SiteID=STAID) %>%
    filter(!is.na(chnm),
           CAS != "---" | is.na(CAS),
           CAS != "-" | is.na(CAS))
  
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

clean_names <- function(cas_df){
  
  ignore_totals <- c(#"Total Pcbs",
                     "Total Pcbs in Mg/l",
                     "Total Oc Pesticides",
                     "Tcpp_isomer","Tcpp Isomer")
  
  cas_final =  cas_df %>%
    mutate(chnm = tools::toTitleCase(chnm)) %>% 
    filter(!(chnm %in% ignore_totals))
  
  cas_final$CAS[cas_final$chnm == "Total Pcbs"] <- "1336-36-3"
  cas_final$chnm[cas_final$CAS == "1336-36-3"] <- "Total PCBs"
  
  cas_final$CAS[cas_final$chnm == "Bupropion"] <- "34841-39-9"
  cas_final$CAS[cas_final$chnm == "Buproprion"] <- "34841-39-9" #Spelled wrong in raw data
  cas_final$chnm[cas_final$CAS == "34841-39-9"] <- "Bupropion"

  # cas_final$chnm[cas_final$CAS == "34911-55-2"] <- "Bupropion hydrochloride"
  
  
  cas_final$CAS[cas_final$chnm == "Nadolol"] <- "42200-33-9"
  cas_final$CAS[cas_final$chnm == "Omeprazole + Esomprazole"] <- "73590-58-6"
  cas_final$CAS[cas_final$chnm == "Tris(1,3-Dichloro-2-Propyl)Phosphate (t"] <- "13674-87-8"
  
  cas_final$CAS[cas_final$chnm == "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"] <- "13674-87-8"
  cas_final$CAS[cas_final$CAS == "26248-87-3"] <- "13674-84-5" #2 versions of TDCPP
  cas_final$CAS[cas_final$CAS == "51805-45-9"] <- "115-96-8" # 2 versions of TCEP

  cas_final$CAS[cas_final$chnm == "Tri(2-chloroethyl) phosphate (TCEP)"] <- "115-96-8"

  cas_final$CAS[cas_final$CAS == "30306-93-5"] <- "77-93-0"
  cas_final$CAS[cas_final$CAS == "53-07-3"] <- "53-70-3"
  
  cas_final$chnm[cas_final$CAS == "77-93-0"] <- "Ethyl Citrate"

  cas_final$chnm[cas_final$CAS == "73590-58-6"] <- "Omeprazole + Esomprazole"
  cas_final$chnm[cas_final$CAS == "101-20-2"] <- "3,4,4'-Trichlorocarbanilide"
  cas_final$chnm[cas_final$CAS == "115-96-8"] <- "Tri(2-chloroethyl)phosphate (TCEP)"
  cas_final$chnm[cas_final$chnm == "Triethyl Phosphate (Tep)"] <- "Triethyl Phosphate (TEP)"
  cas_final$chnm[cas_final$chnm == "Deet"] <- "N,N-diethyltoluamide (DEET)"
  cas_final$chnm[cas_final$chnm == "N,n-Diethyltoluamide (Deet)"] <- "N,N-diethyltoluamide (DEET)"
  cas_final$chnm[cas_final$chnm =="Diethylhexylphthalate (Dehp)"] <- "Diethylhexylphthalate (DEHP)"
  cas_final$chnm[cas_final$chnm == "Tris(2-Ethylhexyl)Phosphate (Tehp)"] <- "Tris(2-ethylhexyl)phosphate (TEHP)"
  cas_final$chnm[cas_final$chnm == "Triphenyl Phosphate (Tpp)"] <- "Triphenyl Phosphate (TPP)"

  cas_final$chnm[cas_final$chnm == "Tris(1-Chloro-2-Propyl)Phosphate (Tcpp)"] <- "Tris(1-chloro-2-propyl)phosphate (TCPP)"
  cas_final$chnm[cas_final$chnm == "Tris(2-Butoxyethyl)Phosphate (Tbep)"] <- "Tris(2-butoxyethyl)phosphate (TBEP)"
  cas_final$chnm[cas_final$chnm == "Tcpp"] <- "Tris(1-chloro-2-propyl)phosphate (TCPP)"
  cas_final$chnm[cas_final$chnm == "Tbep"] <- "Tris(2-butoxyethyl)phosphate (TBEP)"
  cas_final$chnm[cas_final$chnm == "Tdcpp"] <- "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"
  cas_final$chnm[cas_final$chnm == "Tris(1,3-Dichloro-2-Propyl)Phosphate (t"] <- "Tris(1,3-dichloro-2-propyl)phosphate (TDCPP)"
  cas_final$chnm[cas_final$chnm == "Total Pcbs"] <- "Total PCBs"
  cas_final$chnm[cas_final$chnm == "O,p'-Ddd"] <- "o,p'-DDD"
  cas_final$chnm[cas_final$chnm == "P,p'-Ddd"] <- "p,p'-DDD"
  cas_final$chnm[cas_final$chnm == "Pentachloroanisole (Pca)"] <- "Pentachloroanisole (PCA)"
  cas_final$chnm[cas_final$chnm == "Tributyl Phosphate (Tbp)"] <- "Tributyl Phosphate (TBP)"
  cas_final$chnm[cas_final$chnm == "Hydrochlorothiazide (Hctz)"] <- "Hydrochlorothiazide (HCTZ)"
  cas_final$chnm[cas_final$chnm == "O,p'-Ddt"] <- "o,p'-DDT"
  cas_final$chnm[cas_final$chnm == "O,p'-Ddt"] <- "o,p'-DDT"
  cas_final$chnm[cas_final$chnm == "P,p'-Dde"] <- "p,p'-DDE"
  cas_final$chnm[cas_final$chnm == "P,p'-Ddt"] <- "p,p'-DDT"
  cas_final$chnm[cas_final$chnm == "O,p'-Dde"] <- "o,p'-DDE"
  cas_final$chnm[cas_final$chnm == "Indeno[1,2,3-Cd]pyrene"] <- "Indeno[1,2,3-cd]pyrene"
  cas_final$chnm[cas_final$chnm == "Benzo(a)Pyrene"] <- "Benzo(a)pyrene"
  cas_final$chnm[cas_final$chnm == "beta-Bhc"] <- "beta-Hexachlorocyclohexane"
  cas_final$chnm[cas_final$chnm == "P,p'-Methoxychlor"] <- "p,p'-Methoxychlor"
  cas_final$chnm[cas_final$chnm == "alpha-Bhc"] <- "alpha-Hexachlorocyclohexane"
  cas_final$chnm[cas_final$chnm == "Benzo[b]naphtho[2,1-D]thiophene"] <- "Benzo[b]naphtho[2,1-d]thiophene"
  cas_final$chnm[cas_final$chnm == "Dibenzo[a,h]anthracene"] <- "Dibenz[a,h]anthracene"
  cas_final$chnm[cas_final$chnm == "p-Tert-Octylphenol"] <- "p-tert-Octylphenol"
  cas_final$chnm[cas_final$CAS =="26248-87-3"] <- "Tri(chloropropyl) phosphate"
  cas_final$chnm[cas_final$chnm == "Hexachlorobenzene (Hcb)"] <- "Hexachlorobenzene (HCB)"
  
  cas_final$chnm[cas_final$chnm == "Celestolide (Adbi)"] <- "Celestolide (ADBI)"
  cas_final$chnm[cas_final$chnm == "Traseolide (Atii)"] <- "Traseolide (ATII)"
  cas_final$chnm[cas_final$chnm == "Phantolide (Ahmi)"] <- "Phantolide (AHMI)"
  cas_final$chnm[cas_final$chnm == "Galaxolide (Hhcb)"] <- "Galaxolide (HHCB)"
  cas_final$chnm[cas_final$chnm == "Tonalide (Ahtn)"] <- "Tonalide (AHTN)"
  cas_final$chnm[cas_final$chnm == "Para-Cresol"] <- "para-Cresol"
  cas_final$chnm[cas_final$chnm == "Endosulfan-Ii"]  <- "Endosulfan-II"
  cas_final$chnm[cas_final$CAS == "77-93-0"] <- "Triethyl Citrate"
  
  cas_final$chnm[cas_final$CAS == "101-20-2"] <- "3,4,4'-Trichlorocarbanilide"
  
  cas_final$chnm[grep("Cis-", cas_final$chnm)] <- gsub(pattern = "Cis-",
                                                       replacement = "cis-",
                                                       cas_final$chnm[grep("Cis-", cas_final$chnm)])
  cas_final$chnm[grep("Trans-", cas_final$chnm)] <- gsub(pattern = "Trans-",
                                                         replacement = "trans-",
                                                         cas_final$chnm[grep("Trans-", cas_final$chnm)])
  
  cas_final$chnm[grep("Pbde-", cas_final$chnm)] <- gsub(pattern = "Pbde-",
                                                        replacement = "PBDE-",
                                                        cas_final$chnm[grep("Pbde-", cas_final$chnm)])
  cas_final$chnm[grep("-Bhc", cas_final$chnm)] <- gsub(pattern = "-Bhc",
                                                        replacement = "-BHC",
                                                        cas_final$chnm[grep("-Bhc", cas_final$chnm)])
  cas_final$chnm[grep("Alpha-", cas_final$chnm)] <- gsub(pattern = "Alpha-",
                                                       replacement = "alpha-",
                                                       cas_final$chnm[grep("Alpha-", cas_final$chnm)])
  cas_final$chnm[grep("Beta-", cas_final$chnm)] <- gsub(pattern = "Beta-",
                                                         replacement = "beta-",
                                                         cas_final$chnm[grep("Beta-", cas_final$chnm)])
  cas_final$chnm[grep("Delta-", cas_final$chnm)] <- gsub(pattern = "Delta-",
                                                        replacement = "delta-",
                                                        cas_final$chnm[grep("Delta-", cas_final$chnm)])
  

  
  return(cas_final)
  
}

clean_cas <- function(cas_df){
  
  cas_final =  cas_df %>%
    filter(!duplicated(CAS)) 
  
  cas_final <- clean_names(cas_final)

  return(cas_final)
}