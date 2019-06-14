loadd(chemicalSummary)
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

chemicalSummaryPriority <- filter(chemicalSummary, endPoint %in% priority_endpoints)

graphData <- chemicalSummaryPriority %>%
  group_by(site,date,endPoint) %>%
  summarise(sumEAR=sum(EAR)) %>%
  group_by(endPoint) %>%
  summarise(medianEAR = median(sumEAR),
            maxEAR = max(sumEAR)) 

chem_ns <- chemicalSummaryPriority %>%
  left_join(AOP, by="endPoint") %>%
  group_by(endPoint) %>%
  summarize(chems = list(unique(as.character(chnm))),
            AOPs = list(unique(ID)),
            sites = list(unique(as.character(shortName))))


aop_sum <- graphData %>%
  left_join(chem_ns, by="endPoint") %>%
  arrange(desc(maxEAR))

dir.create("data/supplemental", showWarnings = FALSE)

data.table::fwrite(aop_sum, file = "data/supplemental/AOP_priority.txt", sep = "\t")
