library(shiny)
library(toxEval)
library(tidyverse)
library(here)
source(here("R/mixtures/mix_script.R"))
source(here("R/plot_endpoints_facet.R"))

drake::loadd(chemicalSummary)
site_info <- drake::readd(site_info)
chem_info <- drake::readd(chem_info) %>% 
    rename(Chemical = chnm)
chem_data <- drake::readd(all_data_fixed_cas)

end_point_info <- end_point_info

ep_summary <- chemicalSummary %>% 
    mutate(chnm = as.character(chnm)) %>% 
    group_by(chnm, CAS) %>% 
    summarise(n_endPoints = length(unique(endPoint)),
              n_sites_det = length(unique(shortName[EAR > 0])))

# What we lose:
no_ep_chem_info <- chem_info %>% 
    filter(!(CAS %in% unique(ep_summary$CAS)))

no_ep_chem_data = chem_data %>% 
    filter(CAS %in% no_ep_chem_info$CAS)

tox_list_no_ep <- list(chem_data = no_ep_chem_data,
                       chem_info = no_ep_chem_info,
                       chem_site = site_info) %>% 
    as.toxEval()

# Chemicals we lose with toxCast:
chemicalSummary_no_ep <- get_concentration_summary(tox_list = tox_list_no_ep)

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

missing_endpoints_aop <- chemicalSummary %>% 
    left_join(select(AOP_crosswalk, endPoint, ID), by="endPoint") %>% 
    filter(is.na(ID)) %>% 
    select(endPoint) %>% 
    distinct() %>% 
    pull(endPoint)

chemicalSummary_aop_split <- chemicalSummary

chemicalSummary_aop_split$guide_side <- "Keep"
chemicalSummary_aop_split$guide_side[chemicalSummary_aop_split$endPoint %in% missing_endpoints_aop] <- "Lose"

chemicalSummary_aop_split <- chemicalSummary_aop_split %>% 
    bind_rows(mutate(chemicalSummary, guide_side = "All")) %>% 
    mutate(guide_side = factor(guide_side, levels = c("All","Keep","Lose")))

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

missing_endpoints_genes <- chemicalSummary %>% 
    left_join(gene, by="endPoint") %>% 
    filter(is.na(gene)) %>% 
    select(endPoint) %>% 
    distinct() %>% pull(endPoint)

gene_summary <- chemicalSummary %>% 
    left_join(gene, by="endPoint") %>%
    filter(!is.na(gene)) %>% 
    select(-endPoint) %>% 
    mutate(endPoint = gene)

chemicalSummary_gene_split <- chemicalSummary

chemicalSummary_gene_split$guide_side <- "Keep"
chemicalSummary_gene_split$guide_side[chemicalSummary_gene_split$endPoint %in% missing_endpoints_genes] <- "Lose"

chemicalSummary_gene_split <- chemicalSummary_gene_split %>% 
    bind_rows(mutate(chemicalSummary, guide_side = "All")) %>% 
    mutate(guide_side = factor(guide_side, levels = c("All","Keep","Lose")))

#######################################
# Panther
panther <- data.table::fread(here("panther_data/joined_genes.csv")) %>% 
    data.frame() %>% 
    select(gene = gene_abbr,
           pathway_accession,
           pathway_name) %>% 
    left_join(gene,  by="gene")

missing_ep_panther <- chemicalSummary_gene_split %>% 
    filter(guide_side == "Keep") %>% 
    select(-guide_side) %>% 
    left_join(panther, by="endPoint") %>%  
    filter(pathway_accession == "") %>% 
    select(endPoint) %>% distinct() %>% pull(endPoint)

panther_summary <- chemicalSummary %>% 
    left_join(panther, by="endPoint") %>%
    filter(!is.na(gene)) %>% 
    filter(pathway_accession != "") %>% 
    select(-endPoint) %>% 
    mutate(endPoint = pathway_name)

chemicalSummary_panther_split <- chemicalSummary_gene_split %>% 
    filter(guide_side == "Keep") %>% 
    select(-guide_side)

chemicalSummary_panther_split$guide_side <- "Keep"
chemicalSummary_panther_split$guide_side[chemicalSummary_panther_split$endPoint %in% missing_ep_panther] <- "Lose"

chemicalSummary_panther_split <- chemicalSummary_panther_split %>% 
    bind_rows(chemicalSummary_gene_split %>% 
                   filter(guide_side == "Keep") %>% 
                   mutate(guide_side = "All from Genes")) %>% 
    mutate(guide_side = factor(guide_side, levels = c("All from Genes","Keep","Lose")))

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

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    output$tox_missing <- renderPlot({ 
        conc_plot <- plot_tox_boxplots(chemicalSummary_no_ep,
                                       category = "Chemical", 
                                       hit_threshold = input$ear_thresh,
                                       title = "Chemicals lost by ToxCast",
                                       x_label = "Concentration [ug/L]")
        return(conc_plot)
    })
    
    output$aop_missing <- renderPlot({ 
        plot_chems <- input$chemical == "Chemicals"
        if(plot_chems){
            aop_plot <- plot_chemical_boxplots(chemicalSummary_aop_split,
                                               guide_side, 
                                               hit_threshold = input$ear_thresh, 
                                               plot_ND = FALSE) +
                ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
        } else {
            
            ch <- chemicalSummary_aop_split %>% 
                filter(guide_side %in% c("Lose","Keep")) 
            
            aop_plot <- plot_tox_endpoints2(ch, guide_side, 
                                             hit_threshold = input$ear_thresh,
                                             title = "Lost EARs going to AOP",
                                             category = "Chemical") +
                ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
        }

        return(aop_plot)
    })
    
    output$gene_missing <- renderPlot({ 
        plot_chems <- input$chemicalGene == "Chemicals"
        
        if(plot_chems){
            gene_plot <- plot_chemical_boxplots(chemicalSummary_gene_split,
                                                guide_side, 
                                                hit_threshold = input$ear_thresh, 
                                                plot_ND = FALSE) +
                ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
        } else {
            
            
            ch <- chemicalSummary_gene_split %>% 
                filter(guide_side %in% c("Lose","Keep")) 
            
            gene_plot <- plot_tox_endpoints2(ch, guide_side, 
                                                hit_threshold = input$ear_thresh,
                                                title = "Lost EARs going to genes",
                                                category = "Chemical") +
                ggplot2::facet_grid(. ~ guide_side, scales = "free_x")

        }
        

        return(gene_plot)
    })
    
    output$panther_missing <- renderPlot({ 

        if(input$chemicalPath == "Chemicals"){
            panther_plot <- plot_chemical_boxplots(chemicalSummary_panther_split,
                                                   guide_side, 
                                                   hit_threshold = input$ear_thresh,
                                                   plot_ND = FALSE) +
                ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
        } else if (input$chemicalPath == "Endpoints") {
            
            ch <- chemicalSummary_panther_split %>% 
                filter(guide_side %in% c("Lose","Keep")) 
            
            panther_plot <- plot_tox_endpoints2(ch, guide_side, 
                                            hit_threshold = input$ear_thresh,
                                            title = "Lost EARs going to Panther",
                                            category = "Chemical") +
                ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
            
        } else {
            
            ch <- chemicalSummary_gene_split %>% 
                filter(guide_side == "Keep") %>% 
                select(-guide_side) %>% 
                left_join(gene, by="endPoint") %>% 
                left_join(select(panther, -gene), by="endPoint") %>%
                select(-endPoint) %>% 
                rename(endPoint = gene) 
            
            ch$guide_side <- "Keep"
            ch$guide_side[ch$pathway_accession == ""] <- "Lose"

            panther_plot <- plot_tox_endpoints2(ch, guide_side, 
                                                hit_threshold = input$ear_thresh,
                                                title = "Lost Genes going to Panther",
                                                category = "Chemical") +
                ggplot2::facet_grid(. ~ guide_side, scales = "free_x")
                
        }
        return(panther_plot)
    })
    
    overall_mixes <- reactive({
        ear_thresh <- input$ear_thresh
        EARsum_endpoint <- sum_endpoints(chemicalSummary,
                                         ear_cutoff = ear_thresh)
        top_mixes <- all_mixes_fn(EARsum_endpoint, ear_thresh)
    })
    
    output$overallSummary <- DT::renderDataTable({
        
        n_site_thresh <- input$n_sites
        
        top_mixes <- overall_mixes() %>% 
            filter(n_sites >= n_site_thresh)
        
        df <- join_everything  %>% 
            filter(endPoint %in% top_mixes$endPoint)
        
        return(df)
    })
    
    ep_mix_all <- reactive({
        ear_thresh <- input$ear_thresh
        EARsum_endpoint <- sum_endpoints(chemicalSummary,
                                         ear_cutoff = ear_thresh)
        top_mixes <- all_mixes_fn(EARsum_endpoint, ear_thresh)
    })
    
    output$epMixes <- DT::renderDataTable({ 
        
        n_site_thresh <- input$n_sites
        
        top_mixes <- ep_mix_all() %>% 
            filter(n_sites >= n_site_thresh) %>% 
            arrange(desc(n_sites))
        
        top_mixes_dt <- DT::datatable(top_mixes, caption = "ToxCast",
                                      rownames = FALSE,
                                      options = list(scrollX = TRUE,
                                                     pageLength = 5))
        
        return(top_mixes_dt)
    })
    
    
    aop_mix_all <- reactive({
        ear_thresh_aop <- input$ear_thresh
        EARsum_endpoint_aop <- sum_endpoints(aop_summary,
                                             ear_cutoff = ear_thresh_aop)
        top_mixes <- all_mixes_fn(EARsum_endpoint_aop,
                                  ear_thresh_aop)
    })
    
    
    output$aopMixes <- DT::renderDataTable({
      
        n_site_thresh_aop <- input$n_sites

        top_mixes_aop <- aop_mix_all() %>% 
            filter(n_sites >= n_site_thresh_aop) %>% 
            rename(`AOP ID` = endPoint) %>% 
            arrange(desc(n_sites))
        
        top_mixes_aop_dt <- DT::datatable(top_mixes_aop, caption = "AOP",
                                          rownames = FALSE,
                                          options = list(scrollX = TRUE,
                                                         pageLength = 5)) 
        return(top_mixes_aop_dt)
    })
    
    gene_mixes <- reactive({
        ear_thresh_gene <- input$ear_thresh
        EARsum_endpoint_gene <- sum_endpoints(gene_summary,
                                              ear_cutoff = ear_thresh_gene)
        gene_mixes <- all_mixes_fn(EARsum_endpoint_gene,
                                   ear_thresh_gene)
        
    })
    
    output$geneMixes <- DT::renderDataTable({
        
        n_site_thresh_gene <- input$n_sites
        
        top_mixes_gene <- gene_mixes() %>%
            filter(n_sites > n_site_thresh_gene) %>% 
            rename(Gene = endPoint) %>% 
            arrange(desc(n_sites))
        
        top_mixes_gene_dt <- DT::datatable(top_mixes_gene, caption = "Gene",
                                              rownames = FALSE,
                                              options = list(scrollX = TRUE,
                                                             pageLength = 5))
        
        return(top_mixes_gene_dt)
    })
    
    top_panters <- reactive({
        ear_thresh_panther <- input$ear_thresh
        EARsum_endpoint_panther <- sum_endpoints(panther_summary,
                                                 ear_cutoff = ear_thresh_panther)
        top_panters <- all_mixes_fn(EARsum_endpoint_panther,
                     ear_thresh_panther)
        
    })
    
    output$pantherMixes <- DT::renderDataTable({
        
        n_site_thresh_panther <- input$n_sites

        top_mixes_panther <- top_panters() %>% 
                filter(n_sites >= n_site_thresh_panther)
                rename(Panther_Name = endPoint) %>% 
                arrange(desc(n_sites))
        
        top_mixes_panther_dt <- DT::datatable(top_mixes_panther, caption = "Panther",
                                              rownames = FALSE,
                                              options = list(scrollX = TRUE,
                                                             pageLength = 5))
        
        return(top_mixes_panther_dt)
    })
})


