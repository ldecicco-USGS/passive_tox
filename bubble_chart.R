library(toxEval)
library(ggplot2)
library(dplyr)

path_to_file <- 'cleanedData/passive.xlsx' 
tox_list <- create_toxEval(path_to_file)
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

raw_data <- tox_list$chem_data

summary_data <- raw_data %>%
  left_join(select(tox_list$chem_info, CAS, Class), by="CAS") %>%
  group_by(Class, SiteID) %>%
  summarise(totals = length(unique(`Sample Date`)),
            sum_conc = sum(Value,na.rm = TRUE)/length(unique(`Sample Date`)),
            freq = sum(is.na(comment) | comment != "<",na.rm = TRUE)/n()) %>%
  ungroup() %>%
  mutate(what = "Conc",
         rel_sum = sum_conc/max(sum_conc))

summary_EAR <- chemical_summary %>%
  rename(SiteID=site) %>%
  group_by(Class, SiteID) %>%
  summarise(sum_conc = sum(EAR, na.rm = TRUE)/n(),
            freq = sum(EAR != 0, na.rm = TRUE)/n()) %>%
  ungroup() %>%
  mutate(what = "EAR",
         rel_sum = sum_conc/max(sum_conc))

summary_EAR$rel_sum[summary_EAR$sum_conc == 0] <- 0

sum_tots <- bind_rows(summary_data, summary_EAR)

my_breaks <- c(10^-8, 10^-6, 10^-4, 10^-2, 10^0)

bubbles <- ggplot(data = filter(sum_tots, !(Class %in% c("Sterol","PCBs","Fuel","Dye/Pigment",
                                              "Human Drug, Non-prescription")))) +
  geom_point(aes(x = SiteID, y = Class, size = freq, color = rel_sum)) +
  facet_grid(. ~ what, scales = "free") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.y = element_text(color = "black", vjust = 0.2), 
        axis.text.x = element_text(color = "black", vjust = 0, angle = 90),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  scale_size(range=c(2, 6)) +
  scale_color_gradient(low = "white", high = "red", 
                       trans = "log", breaks = my_breaks)

ggsave(bubbles, filename = "plots/bubbles.pdf", height = 9, width = 11)
