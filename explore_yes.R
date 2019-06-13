library(drake)
library(tidyverse)
library(googledrive)
library(readxl)
library(data.table)
library(toxEval)
library(openxlsx)

source("R/analyze/data_reader.R")
loadd(cas_df)

YES_2014 = generic_file_opener("data/raw/YES_2014.xlsx",
                               cas_df,
                               n_max = 2, 
                               sheet = "YES",
                               site_sheet = "site info",
                               year = 2014)
YES_2010 = generic_file_opener("data/raw/YES_2010.xlsx",
                               cas_df,
                               n_max = 2, 
                               sheet = "YES",
                               site_sheet = "site info",
                               year = 2010,
                               skip_site = 2)

yes_data <- bind_rows(YES_2010,
                      YES_2014)
yes_data$CAS <- 1
yes_info <- data.frame(CAS = 1,
                       Class = "YES",
                       chnm = "YES",
                       stringsAsFactors = FALSE)

loadd(sites)
loadd(chemicalSummary)

eps <- c("ATG_ERE_CIS_up","ATG_Era_Trans_up","TOX21_Era_BLA_Agonist_ch1",
         "TOX21_Era_BLA_Agonist_ch2","TOX21_Era_LUC_BG1_Agonis")

chem_sum_eps <- dplyr::filter(chemicalSummary, endPoint %in% eps)

eps_data <- chem_sum_eps %>%
  group_by(site, date, endPoint) %>%
  summarise(Value = sum(EAR, na.rm = TRUE)) %>%
  ungroup()

yes_data <- yes_data %>%
  select(site = SiteID, date=`Sample Date`, YES=Value) 

all_data <- eps_data %>%
  left_join(yes_data, by=c("site","date")) %>%
  filter(!is.na(YES))

library(ggplot2)
library(quantreg)


yes_plot <- ggplot(all_data, aes(x=Value, y=YES)) +
  geom_point() +
  facet_grid(. ~ endPoint) + 
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_quantile(quantiles = c(0.1, 0.5,0.9))

yes_plot

ggsave(yes_plot, filename = "plots/yes.png")


new_order <- chem_sum_eps %>%
  mutate(chnm = as.character(chnm)) %>%
  group_by(chnm) %>%
  summarize(med = median(EAR[EAR!=0], na.rm = TRUE)) %>%
  arrange(med) 

if(any(is.na(new_order$med))){
  orderedLevels <- c(new_order$chnm[is.na(new_order$med)],
                     new_order$chnm[!is.na(new_order$med)])
} else {
  orderedLevels <- new_order$chnm
}

chem_sum_eps$chnm <- factor(as.character(chem_sum_eps$chnm), 
                            levels = orderedLevels)

n_fun <- function(x){
  return(data.frame(y = -9,
                    label = length(x)))
}

chem_contrib <- ggplot(chem_sum_eps, aes(x=chnm,  y=EAR)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
  coord_flip() +
  xlab("") +
  ggtitle("EAR by chemical for ATG_ERE_CIS_up")

chem_contrib
ggsave(chem_contrib, filename = "plots/yes_chems.png", width = 8, height = 6)

chem_yes_data <- chem_sum_eps %>%
  left_join(yes_data, by=c("site","date")) %>%
  filter(!is.na(YES))

chem_yes_data_filtered <- chem_yes_data %>%
  filter(EAR > 0, YES > 0)

chem_facets <- ggplot(data = chem_yes_data_filtered) +
  geom_point(aes(x=EAR,y=YES)) +
  facet_wrap(. ~ chnm) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")
  
ggsave(chem_facets, filename = "plots/yes_chems_scatter.png", width = 11, height = 9)



