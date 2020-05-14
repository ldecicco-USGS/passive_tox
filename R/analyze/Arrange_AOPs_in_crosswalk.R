#Find AOPs
library(tidyverse)

AOPs <- read.csv(file.path("data","raw","AOP_crosswalk.csv"))
AOPs <- AOPs %>% mutate(row_num = row_number())
AOPs <- AOPs[,c(8,1:7)]
AOPs_in_included <- c(7,
                     8,
                     11,
                     18,
                     21,
                     25,
                     41,
                     57,
                     60,
                     122,
                     123,
                     131,
                     150,
                     153)
AOPs_in_passive <- AOPs[which(AOPs$AOP.. %in% AOPs_in_included),] %>%
  arrange(AOP..)
