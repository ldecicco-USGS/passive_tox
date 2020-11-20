library(tidyverse)
library(toxEval)
library(openxlsx)

plot_mixes <- function(){
  path_to_data <- Sys.getenv("PASSIVE_PATH")
  
  options(dplyr.summarise.inform = FALSE)
  
  # shiny::runApp("apps/Mixture_Exploration/")
  
  # Data setup is done in drake plan
  # It's better to run separately, and shouldn't
  # be needed to run this script
  
  source(file = "read_chemicalSummary.R")
  
  unique(tox_list$chem_data$CAS)[which(!(unique(tox_list$chem_data$CAS)) %in% tox_list$chem_info$CAS)]
  
  tox_list$exclusions <- tox_list$exclusions %>% 
    filter(!is.na(CAS) & !is.na(endPoint))
  
  ###################################################
  # Mixtures stuff:
  EAR_thresh <- 0.001
  site_thresh_percent <- 10
  n_sites <- 7 
  TQ_thresh <- 0.1
  
  source(file = "R/mixtures/mix_script.R")
  source(file = "R/mixtures/prepare_mixture_data.R")
  
  mix_df <- get_final_mixtures(chemicalSummary,
                               EAR_thresh,
                               site_thresh_percent, tox_list)
  
  mix_slim <- mix_df %>% 
    select(Genes, siteNames) %>% 
    distinct() %>% 
    mutate(siteNames = strsplit(siteNames, "\\|"))
  
  mix_sites <- tox_list$chem_site %>% 
    select(`Short Name`, site_grouping)
  
  for(i in 1:nrow(mix_slim)){
    mix_sites[[mix_slim$Genes[i]]] <- mix_sites$`Short Name` %in% unlist(mix_slim$siteNames[i])
  }
  
  # Table was corrected by hand to do 1 row per gene:
  mix_sites$NR1I2 <- mix_sites$`Short Name` %in% unique(unlist(mix_slim$siteNames[4:6]))
  mix_sites$CYP2B6 <- mix_sites$`Short Name` %in% unique(unlist(mix_slim$siteNames[10:11]))
  
  mix_sites_long <- mix_sites %>% 
    pivot_longer(cols = c(-`Short Name`, -site_grouping),
                 names_to = "Gene Symbols")
  
  sites_ordered <- c("StLouis","Nemadji","WhiteWI","Bad",
                     "Montreal","PresqueIsle","Pigeon","Ontonagon",
                     "Sturgeon","Tahquamenon",
                     "Manistique","Escanaba","Ford","Menominee",
                     "Peshtigo","Oconto","Fox","Manitowoc",
                     "Sheboygan #4","Sheboygan #3","Sheboygan #2","Sheboygan #1",
                     "MilwaukeeMouth","IndianaHC #2", "IndianaHC #1",
                     "Burns","StJoseph","PawPaw","Kalamazoo",
                     "GrandMI #1","GrandMI #2","GrandMI #3","GrandMI #4",
                     "Muskegon","WhiteMI","PereMarquette","Manistee",
                     "Indian","Cheboygan","ThunderBay","AuSable",
                     "Rifle","Saginaw","BlackMI","Clinton",
                     "Rouge","HuronMI","Raisin",
                     "Maumee","Portage","Sandusky","HuronOH",
                     "Vermilion","BlackOH","Rocky","Cuyahoga",
                     "GrandOH","Ashtabula","Cattaraugus","Buffalo",
                     "Tonawanda","Genesee","GeneseeDock","Oswego","BlackNY",
                     "Oswegatchie","Grass","Raquette","StRegis")
  
  mix_sites_long$`Short Name` <- factor(mix_sites_long$`Short Name`, levels = rev(sites_ordered))
  mix_sites_long$site_grouping <- factor(mix_sites_long$site_grouping, levels = c("Lake Superior",
                                                                                  "Lake Michigan",
                                                                                  "Lake Huron",
                                                                                  "Lake Erie",
                                                                                  "Lake Ontario"))
  
  plot_out <- ggplot(data = mix_sites_long) +
    geom_tile(aes(x = `Gene Symbols`, y = `Short Name`, fill = value)) +
    facet_grid(site_grouping ~ ., scales = "free", space = "free") +
    theme_minimal() +
    scale_fill_manual("Detected mixture", labels = c("No", "Yes"),
                      values = c("transparent", "steelblue")) +
    theme(legend.position=c(.9,.93),
          axis.title.y = element_blank())
  
  return(plot_out)
}
