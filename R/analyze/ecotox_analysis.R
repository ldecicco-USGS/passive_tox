tox <- read.delim(file = file.path(Sys.getenv("PASSIVE_PATH"),"ECOTOX","DEHP.txt"),sep="|", stringsAsFactors = FALSE)

tox$X.Conc.1.Mean..Standardized.. <- as.numeric(tox$X.Conc.1.Mean..Standardized..)
tox$X.Conc.Min.1..Standardized.. <- as.numeric(tox$X.Conc.Min.1..Standardized..)
tox.test <- mutate(endpoint = min(tox[,c("X.Conc.1.Mean..Standardized..","X.Conc.Min.1..Standardized..")]))

tox_fw <- tox %>%
  filter(Media.Type == "Fresh water",
         Conc.1.Type..Standardized.. == "Active ingredient")

boxplot(tox_fw$X.Conc.1.Mean..Standardized..)


        
        names(tox)
        