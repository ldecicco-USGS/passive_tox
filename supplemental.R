library(drake)
library(toxEval)
library(tidyverse)

loadd(chemicalSummary)

length(unique(chemicalSummary$endPoint))

end_point_info <- end_point_info

ep_we_use <- chemicalSummary %>% 
  select(endPoint) %>% 
  distinct() %>% 
  left_join(select(end_point_info, 
                   endPoint = assay_component_endpoint_name,
                   assay = assay_source_long_name), by="endPoint")

any(ep_we_use$assay == "NHEERL Padilla Lab")

zebra <- chemicalSummary %>% 
  left_join(ep_we_use, by="endPoint") %>% 
  filter(assay == "NHEERL Padilla Lab")

######################################
# AOP Plot details
loadd(AOP)


threshold = 0.001
siteThreshold = 10

endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= threshold) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites)) %>%
  filter(numSites >= siteThreshold)

priority_endpoints <- endpoints_sites_hits$endPoint

AOP_list <- AOP %>% 
  filter(endPoint %in% priority_endpoints) %>% 
  group_by(endPoint) %>% 
  summarize(n_AOPs = length(unique(ID)),
            AOPs = list(unique(ID)))

nChems <- chemicalSummary %>%
  filter(endPoint %in% priority_endpoints) %>% 
  group_by(endPoint) %>%
  summarize(nChems = length(unique(chnm)),
            Chemicals = list(unique(as.character(chnm))))

graphData <- chemicalSummary %>%
  group_by(shortName,date,endPoint) %>%
  summarise(sumEAR=sum(EAR)) %>%
  group_by(shortName, endPoint) %>%
  summarise(maxEAR=max(sumEAR))

orderColsBy <- graphData %>%
  filter(endPoint %in% priority_endpoints) %>% 
  group_by(endPoint) %>%
  summarise(median = quantile(maxEAR[maxEAR != 0],0.5)) %>%
  arrange(median)

orderedLevelsEP <- orderColsBy$endPoint

aop_summary <- graphData %>%
  ungroup() %>% 
  filter(maxEAR > 0) %>%
  filter(endPoint %in% priority_endpoints) %>% 
  group_by(endPoint) %>%
  summarise(medianEAR = median(maxEAR),
            maxEAR = max(maxEAR),
            n_sites = length(unique(shortName)),
            sites = list(unique(shortName))) %>%
  ungroup() %>%
  left_join(AOP_list, by="endPoint") %>%
  left_join(nChems, by="endPoint")

aop_summary$endPoint <- factor(aop_summary$endPoint, levels = orderedLevelsEP)
aop_summary <- arrange(aop_summary, desc(endPoint))

aop_summary$AOPs[sapply(aop_summary$AOPs, is.null)] <- NA

data.table::fwrite(aop_summary, "data/supplemental/aop_summary.txt", sep = "\t")
