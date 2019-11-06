library(shiny)
library(toxEval)
library(tidyverse)
library(here)
source(here("R/mixtures/mix_script.R"))

drake::loadd(chemicalSummary)

ACC_all = get_ACC(drake::readd(chem_info)$CAS) %>%
    remove_flags()

end_point_info <- end_point_info
cleaned_ep <- clean_endPoint_info(end_point_info)

filtered_ep <- filter_groups(cleaned_ep, 
                            groupCol = 'intended_target_family',
                            assays = c('ATG','NVS','OT',
                                       'TOX21','CEETOX','APR',
                                       'CLD','TANGUAY','NHEERL_PADILLA',
                                       'NCCT_SIMMONS','ACEA'),
                            remove_groups = c('Background Measurement',
                                              'Undefined'))

ep_summary <- chemicalSummary %>% 
    mutate(chnm = as.character(chnm)) %>% 
    group_by(chnm, CAS) %>% 
    summarise(n_endPoints = length(unique(endPoint)),
              n_sites_det = length(unique(shortName[EAR > 0])))

# What we lose:
no_ep_chem_info <- drake::readd(chem_info) %>% 
    filter(!(CAS %in% unique(ep_summary$CAS)))

no_ep_chem_data = drake::readd(all_data_fixed_cas) %>% 
    filter(CAS %in% no_ep_chem_info$CAS)


tox_list_no_ep <- list(chem_data = no_ep_chem_data,
                       chem_info = no_ep_chem_info,
                       chem_site = drake::readd(site_info),
                       benchmarks = drake::readd(benchmarks)) %>% 
    as.toxEval()

chemicalSummary_no_ep <- get_chemical_summary(tox_list = tox_list_no_ep)

#################################################
# AOP
AOP_crosswalk <- data.table::fread("D:/LADData/toxManuscript/AOP_crosswalk_Dec_2018.csv")%>%
    data.frame() %>%
    select(endPoint=Component.Endpoint.Name, ID=AOP.., KE = KE., everything()) %>%
    distinct()

aop_summary <- chemicalSummary %>% 
    left_join(select(AOP_crosswalk, endPoint, ID), by="endPoint") %>% 
    filter(!is.na(ID)) %>% 
    mutate(chnm = as.character(chnm)) %>% 
    group_by(chnm, CAS) %>%
    summarise(n_AOP = length(unique(ID)),
              n_sites_det = length(unique(shortName[EAR > 0])))  

missing_chem <- ep_summary$chnm[!(ep_summary$chnm %in% unique(aop_summary$chnm))]
missing_cas <- ep_summary$CAS[!(ep_summary$CAS %in% unique(aop_summary$CAS))]

# What we lose:
no_aop_chem_info <- drake::readd(chem_info) %>% 
    filter(CAS %in% missing_cas)

no_aop_chem_data = drake::readd(all_data_fixed_cas) %>% 
    filter(CAS %in% no_aop_chem_info$CAS)

tox_list_no_AOPs <- list(chem_data = no_aop_chem_data,
                         chem_info = no_aop_chem_info,
                         chem_site = drake::readd(site_info),
                         exclude = drake::readd(exclude)) %>% 
    as.toxEval()

chemicalSummary_no_aop <- get_chemical_summary(tox_list = tox_list_no_AOPs,
                                               ACC = ACC_all, 
                                               filtered_ep = filtered_ep)
##########################################
# Gene targets
gene_info <- select(end_point_info,
                    endPoint = assay_component_endpoint_name,
                    geneID = intended_target_gene_id,
                    geneName = intended_target_gene_name,
                    geneSymbol = intended_target_gene_symbol)

multiple_genes <- unique(gene_info$geneSymbol[grep(pattern = "\\|",
                                                   x = gene_info$geneSymbol)])
gene_info_wide <- gene_info %>%
    mutate(orig_symbol = geneSymbol) %>% 
    separate(geneSymbol, into = c("a","b","c","d"),
             sep = "\\|") %>% 
    separate(geneID, into = c("IDa","IDb","IDc","IDd"),
             sep = "\\|") %>% 
    separate(geneName, into = c("Namea","Nameb","Namec","Named"),
             sep = "\\|")

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
    mutate(chnm = as.character(chnm)) %>% 
    filter(!is.na(gene)) %>% 
    group_by(chnm, CAS) %>%
    summarise(n_genes = length(unique(gene)),
              n_sites_det = length(unique(shortName[EAR > 0])))

gene_summary <- chemicalSummary %>% 
    left_join(gene, by="endPoint") %>%  
    mutate(chnm = as.character(chnm)) %>% 
    filter(!is.na(gene)) %>% 
    group_by(chnm, CAS) %>%
    summarise(n_genes = length(unique(gene)),
              n_sites_det = length(unique(shortName[EAR > 0])))  

missing_chem_gene <- ep_summary$chnm[!(ep_summary$chnm %in% unique(gene_summary$chnm))]
missing_cas_gene <- ep_summary$CAS[!(ep_summary$CAS %in% unique(gene_summary$CAS))]

no_path_chem_info <- drake::readd(chem_info) %>% 
    filter(CAS %in% missing_cas_gene)

no_path_chem_data = drake::readd(all_data_fixed_cas) %>% 
    filter(CAS %in% no_path_chem_info$CAS)

tox_list_no_path <- list(chem_data = no_path_chem_data,
                         chem_info = no_path_chem_info,
                         chem_site = drake::readd(site_info),
                         exclude = drake::readd(exclude)) %>% 
    as.toxEval()


chemicalSummary_no_gene <- get_chemical_summary(tox_list = tox_list_no_path,
                                                   ACC = ACC_all, 
                                                   filtered_ep = filtered_ep)
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
    group_by(chnm, CAS) %>%
    summarise(n_pathways = length(unique(pathway_accession)),
              n_sites_det = length(unique(shortName[EAR > 0])))  

missing_chem <- ep_summary$chnm[!(ep_summary$chnm %in% unique(panther_summary$chnm))]
missing_cas <- ep_summary$CAS[!(ep_summary$CAS %in% unique(panther_summary$CAS))]

# What we lose:
no_path_chem_info <- drake::readd(chem_info) %>% 
    filter(CAS %in% missing_cas)

no_path_chem_data = drake::readd(all_data_fixed_cas) %>% 
    filter(CAS %in% no_path_chem_info$CAS)

tox_list_no_path <- list(chem_data = no_path_chem_data,
                         chem_info = no_path_chem_info,
                         chem_site = drake::readd(site_info),
                         exclude = drake::readd(exclude)) %>% 
    as.toxEval()


chemicalSummary_no_panther <- get_chemical_summary(tox_list = tox_list_no_path,
                                               ACC = ACC_all, 
                                               filtered_ep = filtered_ep)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    output$tox_missing <- renderPlot({ 
        conc_plot <- plot_tox_boxplots(chemicalSummary_no_ep,
                                       category = "Chemical")
        return(conc_plot)
    })
    
    output$aop_missing <- renderPlot({ 
        aop_plot <- plot_tox_boxplots(chemicalSummary_no_aop,
                                       category = "Chemical")
        return(aop_plot)
    })
    
    output$gene_missing <- renderPlot({ 
        gene_plot <- plot_tox_boxplots(chemicalSummary_no_gene,
                                          category = "Chemical")
        return(gene_plot)
    })
    
    output$panther_missing <- renderPlot({ 
        panther_plot <- plot_tox_boxplots(chemicalSummary_no_panther,
                                      category = "Chemical")
        return(panther_plot)
    })
    
    output$epMixes <- DT::renderDataTable({ 
        ear_thresh <- input$ear_thresh
        n_site_thresh <- input$n_sites
        
        EARsum_endpoint <- sum_endpoints(chemicalSummary,
                                         ear_cutoff = ear_thresh)
        contributing_chems <- calc_contr_chems(EARsum_endpoint)
        top_mixes <- top_mixes_fn(contributing_chems, n_site_thresh)
        
        top_mixes <- top_mixes %>% 
            select(contr_chems_st, endPoint,
                   n_samples, unique_sites)
        
        top_mixes_dt <- DT::datatable(top_mixes, caption = "ToxCast",
                                      rownames = FALSE,
                                      options = list(scrollX = TRUE,
                                                     pageLength = 5))
        
        return(top_mixes_dt)
    })
    
    
    output$aopMixes <- DT::renderDataTable({
        
        ear_thresh_aop <- input$ear_thresh
        n_site_thresh_aop <- input$n_sites
        
        aop_summary <- chemicalSummary %>% 
            left_join(select(AOP_crosswalk, 
                             endPoint, ID), by="endPoint") %>% 
            filter(!is.na(ID)) %>% 
            select(-endPoint) %>% 
            mutate(endPoint = as.character(ID))
        
        EARsum_endpoint_aop <- sum_endpoints(aop_summary,
                                             ear_cutoff = ear_thresh_aop)
        
        contributing_chems_aop <- calc_contr_chems(EARsum_endpoint_aop)
        
        top_mixes_aop <- top_mixes_fn(contributing_chems_aop,
                                      n_site_thresh_aop) %>% 
            select(contr_chems_st, AOP_id = endPoint,
                   n_samples, unique_sites)
        
        top_mixes_aop_dt <- DT::datatable(top_mixes_aop, caption = "AOP",
                                          rownames = FALSE,
                                          options = list(scrollX = TRUE,
                                                         pageLength = 5)) 
        return(top_mixes_aop_dt)
    })
    
    output$geneMixes <- DT::renderDataTable({
        ear_thresh_gene <- input$ear_thresh
        n_site_thresh_gene <- input$n_sites
        
        gene_summary <- chemicalSummary %>% 
            left_join(gene, by="endPoint") %>%
            filter(!is.na(gene)) %>% 
            select(-endPoint) %>% 
            mutate(endPoint = gene)
        
        EARsum_endpoint_gene <- sum_endpoints(gene_summary,
                                              ear_cutoff = ear_thresh_gene)
        
        contributing_chems_gene <- calc_contr_chems(EARsum_endpoint_gene)
        
        top_mixes_gene <- top_mixes_fn(contributing_chems_gene,
                                       n_site_thresh_gene) %>% 
            select(contr_chems_st, Gene_Name = endPoint,
                   n_samples, unique_sites)
        
        top_mixes_gene_dt <- DT::datatable(top_mixes_gene, caption = "Gene",
                                              rownames = FALSE,
                                              options = list(scrollX = TRUE,
                                                             pageLength = 5))
        
        return(top_mixes_gene_dt)
    })
    
    output$pantherMixes <- DT::renderDataTable({
        ear_thresh_panther <- input$ear_thresh
        n_site_thresh_panther <- input$n_sites
        
        panther_summary <- chemicalSummary %>% 
            left_join(panther, by="endPoint") %>%
            filter(!is.na(gene)) %>% 
            filter(pathway_accession != "") %>% 
            select(-endPoint) %>% 
            mutate(endPoint = pathway_name)
        
        EARsum_endpoint_panther <- sum_endpoints(panther_summary,
                                                 ear_cutoff = ear_thresh_panther)
        
        contributing_chems_panther <- calc_contr_chems(EARsum_endpoint_panther)
        
        top_mixes_panther <- top_mixes_fn(contributing_chems_panther,
                                          n_site_thresh_panther) %>% 
            select(contr_chems_st, Panther_Name = endPoint,
                   n_samples, unique_sites)
        
        top_mixes_panther_dt <- DT::datatable(top_mixes_panther, caption = "Panther",
                                              rownames = FALSE,
                                              options = list(scrollX = TRUE,
                                                             pageLength = 5))
        
        return(top_mixes_panther_dt)
    })
})
