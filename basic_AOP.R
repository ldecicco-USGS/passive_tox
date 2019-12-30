library(toxEval)
library(tidyverse)

AOP_info <- readxl::read_xlsx("data/supplemental/SI_6_AOP_relevance With Short AOP name.xlsx", 
                              sheet = "SI_AOP_relevance") %>% 
  group_by(AOP, Relevant, Rationale, `...5`) %>% 
  summarize(`EndPoint(s)` = paste0(`Endpoint(s)`, collapse = ", ")) %>% 
  rename(`Abbreviated AOP description` = `...5`) %>% 
  ungroup()
