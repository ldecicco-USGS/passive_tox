# Combine ECOTOX analysis for toxcast and non toxcast chems

ecotox_toxcast <- read_rds("R/Analyze/Out/ECOTOX_filtered_toxcast.Rds")
ecotox_non_toxcast <- read_rds("R/Analyze/Out/ECOTOX_filtered_non_toxcast.Rds")

keep <- grep("Purity",names(ecotox_toxcast),ignore.case = TRUE,value = TRUE,invert = TRUE)

ecotox <- full_join(ecotox_toxcast[,keep],ecotox_non_toxcast[,keep])

include <- ecotox %>%
  select(Class,CAS,chnm,Species.Scientific.Name., Species.Common.Name, Species.Group, 
         Organism.Lifestage, Organism.Age.Mean.Op., Organism.Age.Mean., Organism.Age.Min.Op.,
         Organism.Age.Min., Organism.Age.Max.Op., Organism.Age.Max., Age.Units, Exposure.Type, 
         Media.Type, Test.Location, Number.of.Doses, 
         Conc.1.Type..Standardized.., Conc.1.Mean.Op..Standardized.., Conc.1.Mean..Standardized.., 
         Conc.1.Min.Op..Standardized.., Conc.Min.1..Standardized.., Conc.1.Max.Op..Standardized.., 
         Conc.1.Max..Standardized.., Conc.1.Units..Standardized.., Conc.2.Type..Standardized.., 
         Conc.2.Mean.Op..Standardized.., Conc.2.Mean..Standardized.., Conc.2.Min.Op..Standardized.., 
         Conc.Min.2..Standardized.., Conc.2.Max.Op..Standardized.., Conc.2.Max..Standardized.., 
         Conc.2.Units..Standardized.., Conc.3.Type..Standardized.., Conc.3.Mean.Op..Standardized.., 
         Conc.3.Mean..Standardized.., Conc.3.Min.Op..Standardized.., Conc.Min.3..Standardized.., 
         Conc.3.Max.Op..Standardized.., Conc.3.Max..Standardized.., Conc.3.Units..Standardized., 
         Effect, Effect.Measurement, Endpoint, Response.Site, 
         Statistical.Significance., Significance.Level.Mean.Op., Significance.Level.Mean., 
         Significance.Level.Min.Op., Significance.Level.Min., Significance.Level.Max.Op., Significance.Level.Max, 
         Observed.Duration.Mean.Op..Days.., Observed.Duration.Mean..Days.., Observed.Duration.Min.Op..Days.., 
         Observed.Duration.Min..Days.., Observed.Duration.Max.Op..Days.., Observed.Duration.Max..Days.., 
         Observed.Duration.Units..Days., Author., Reference.Number., Title., Source., 
         Publication.Year, value)
