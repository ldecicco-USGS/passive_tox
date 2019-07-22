library(drake)
library(tidyverse)
library(googledrive)
library(readxl)
library(data.table)
library(toxEval)
library(openxlsx)

source("R/analyze/data_reader.R")
loadd(cas_df)

eps <- c("ATG_ERE_CIS_up","ATG_ERa_TRANS_up","TOX21_ERa_BLA_Agonist_ch1",
         "TOX21_ERa_BLA_Agonist_ch2","TOX21_ERa_LUC_BG1_Agonist",
         "TOX21_ERa_BLA_Agonist_ratio")
tox_list = create_toxEval("data/clean/passive.xlsx")

ToxCast_ACC <- ToxCast_ACC
x <- ToxCast_ACC %>% 
  filter(endPoint %in% eps,
         CAS %in% tox_list$chem_info$CAS) %>%
  left_join(select(tox_list$chem_info, CAS, chnm), by="CAS") %>%
  left_join(select(end_point_info, endPoint = assay_component_endpoint_name, group = intended_target_family)) %>% 
  filter(group != "background measurement")

x <- x[is.na(x$flags) | 
       !(grepl(pattern = "Borderline active",x = x$flags) |
       grepl(pattern = "Only highest conc",x = x$flags)),]

y <- x %>% 
  arrange(CAS, desc(ACC)) %>% 
  select(-flags, -group) %>% 
  spread(endPoint, ACC)

data.table::fwrite(x, "data/supplemental/yes_endpoints.csv")
data.table::fwrite(y, "data/supplemental/yes_endpoints_wide.csv")

ACC <- get_ACC(tox_list$chem_info$CAS)

all(eps %in% ACC$endPoint)

ACC <- remove_flags(ACC)
all(eps %in% ACC$endPoint)

cleaned_ep <- clean_endPoint_info(end_point_info)
all(eps %in% cleaned_ep$assay_component_endpoint_name)

cleaned_ep$intended_target_family[cleaned_ep$assay_component_endpoint_name %in% c("TOX21_ERa_BLA_Agonist_ch1","TOX21_ERa_BLA_Agonist_ch2")]
cleaned_ep$intended_target_family[cleaned_ep$assay_component_endpoint_name %in% eps]

filtered_ep <- filter_groups(cleaned_ep)
all(eps %in% filtered_ep$endPoint)

eps[!(eps %in% filtered_ep$endPoint)]
eps[(eps %in% filtered_ep$endPoint)]

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

ggsave(yes_plot, filename = "plots/yes.png", width = 10,height = 4)

ave_data <- all_data %>% 
  group_by(site, date) %>% 
  summarise(meanEAR = mean(Value),
            meanYES = mean(YES))

yes_ave_plot <- ggplot(ave_data, aes(x=meanEAR, y=meanYES)) +
  geom_point() +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_quantile(quantiles = c(0.1, 0.5,0.9))

yes_ave_plot

ggsave(yes_ave_plot, filename = "plots/yes_ave.png", width = 6,height = 4)


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
  return(data.frame(y = -10,
                    label = length(x)))
}

chem_contrib <- ggplot() +
  geom_boxplot(data = chem_sum_eps, aes(x=chnm,  y=EAR)) +
  scale_y_continuous(trans = "log10") +
  stat_summary(data = chem_sum_eps, aes(x=chnm,  y=EAR),
               fun.data = n_fun, geom = "text", hjust = 0.5) +
  coord_flip(clip="off") +
  xlab("") +
  theme_bw() +
  facet_grid(. ~ endPoint) 
  # geom_text(aes(y = 1e-10, x = Inf, label = "Hi!"), 
  #           vjust = -0.5) +
  # theme(plot.margin = margin(13, 5.5, 5.5, 5.5, "pt"))

chem_contrib
ggsave(chem_contrib, filename = "plots/yes_chems.png", width = 12, height = 6)

graph_sums <- chem_sum_eps %>%
  group_by(site, chnm, date) %>% 
  summarise(sumEAR = sum(EAR)) %>% 
  group_by(site, chnm) %>% 
  summarise(maxEAR = max(sumEAR))

sumGraph <- ggplot() +
  geom_boxplot(data = graph_sums, aes(x=chnm,  y=maxEAR)) +
  scale_y_continuous(trans = "log10") +
  stat_summary(data = graph_sums, aes(x=chnm,  y=maxEAR),
               fun.data = n_fun, geom = "text", hjust = 0.5) +
  coord_flip(clip="off") +
  xlab("") +
  theme_bw() 

ggsave(sumGraph, filename = "plots/sum_yes_chems.png", width = 12, height = 6)

chem_yes_data <- chem_sum_eps %>%
  left_join(yes_data, by=c("site","date")) %>%
  filter(!is.na(YES))

chem_yes_data_filtered <- chem_yes_data %>%
  filter(EAR > 0, YES > 0)

chem_facets <- ggplot(data = chem_yes_data_filtered) +
  geom_point(aes(x=EAR,y=YES)) +
  facet_wrap(. ~ chnm) +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")
  

ggsave(chem_facets, filename = "plots/yes_chems_scatter.png", width = 11, height = 9)



