library(tidyverse)
library(readxl)
library(toxEval)
library(here)
source(here("R/mixtures/mix_script.R"))

drake::loadd(chemicalSummary)

#################################################
# AOP
AOP_crosswalk <- data.table::fread(here("data/supplemental/AOP_crosswalk_Dec_2018.csv")) %>%
  data.frame() %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP.., KE = KE., everything()) %>%
  distinct()

aop_summary <- chemicalSummary %>% 
  left_join(select(AOP_crosswalk, 
                   endPoint, ID), by="endPoint") %>% 
  filter(!is.na(ID)) %>% 
  select(-endPoint) %>% 
  mutate(endPoint = as.character(ID))


##########################################
# Gene targets
gene_info <- select(end_point_info,
                    endPoint = assay_component_endpoint_name,
                    geneID = intended_target_gene_id,
                    geneName = intended_target_gene_name,
                    geneSymbol = intended_target_gene_symbol)

multiple_genes <- unique(gene_info$geneSymbol[grep(pattern = "\\|",
                                                   x = gene_info$geneSymbol)])

suppressWarnings({
  gene_info_wide <- gene_info %>%
    mutate(orig_symbol = geneSymbol) %>% 
    separate(geneSymbol, into = c("a","b","c","d"),
             sep = "\\|") %>% 
    separate(geneID, into = c("IDa","IDb","IDc","IDd"),
             sep = "\\|") %>% 
    separate(geneName, into = c("Namea","Nameb","Namec","Named"),
             sep = "\\|")  
})


gene_info_long <- gene_info_wide %>% 
  select(endPoint, geneID = IDa, geneSymbol = a, geneName = Namea) %>% 
  rbind(gene_info_wide %>% 
          select(endPoint, geneID = IDb, geneSymbol = b, geneName = Nameb) %>% 
          filter(!is.na(geneID))) %>% 
  rbind(gene_info_wide %>% 
          select(endPoint, geneID = IDc, geneSymbol = c, geneName = Namec) %>% 
          filter(!is.na(geneID))) %>% 
  rbind(gene_info_wide %>% 
          select(endPoint, geneID = IDd, geneSymbol = d, geneName = Named) %>% 
          filter(!is.na(geneID)))

gene <- select(gene_info_long,
               endPoint,
               gene = geneSymbol)

gene_summary <- chemicalSummary %>% 
  left_join(gene, by="endPoint") %>%
  filter(!is.na(gene)) %>% 
  select(-endPoint) %>% 
  mutate(endPoint = gene)


#######################################
# Panther
panther <- data.table::fread(here("panther_data/joined_genes.csv")) %>% 
  data.frame() %>% 
  select(gene = gene_abbr,
         pathway_accession,
         pathway_name) %>% 
  left_join(gene,  by="gene")

panther_summary <- chemicalSummary %>% 
  left_join(panther, by="endPoint") %>%
  filter(!is.na(gene)) %>% 
  filter(pathway_accession != "") %>% 
  select(-endPoint) %>% 
  mutate(endPoint = pathway_name)

join_everything <- chemicalSummary %>% 
  left_join(select(AOP_crosswalk, endPoint, ID), by="endPoint") %>% 
  left_join(gene, by="endPoint") %>% 
  left_join(panther, by=c("endPoint","gene")) %>% 
  select(endPoint, gene, AOP_id = ID, pathway_name) %>% 
  distinct() %>% 
  group_by(endPoint) %>% 
  summarise(AOP_ids = paste(unique(AOP_id[!is.na(AOP_id)]), collapse = ","),
            genes = paste(unique(gene[!is.na(gene)]), collapse = ","),
            pathways = paste(unique(pathway_name[!is.na(pathway_name)]), collapse = ","))

class_key <- chemicalSummary %>% 
  select(chnm, Class) %>% 
  distinct() %>% 
  mutate(chnm = as.character(chnm),
         Class = as.character(Class))
