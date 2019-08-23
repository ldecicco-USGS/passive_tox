graph_chem_data_CAS <- function(chemical_summary, 
                            manual_remove=NULL,
                            mean_logic = FALSE,
                            sum_logic = TRUE){
  
  site <- chnm <- Class <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  chemical_summary <- chemical_summary %>%
    select(-chnm) %>%
    distinct()
  
  if(!sum_logic){
    graphData <- chemical_summary %>%
      dplyr::group_by(site, CAS, Class) %>%
      dplyr::summarise(meanEAR=ifelse(mean_logic,mean(EAR),max(EAR))) %>%
      dplyr::ungroup()     
  } else {
    #With new dplyr...will need to filter out na's in meanEAR
    graphData <- chemical_summary %>%
      dplyr::group_by(site,date,CAS,Class) %>%
      dplyr::summarise(sumEAR=sum(EAR,na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(site, CAS,Class) %>%
      dplyr::summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
      dplyr::ungroup() 
  }
  
  if(!is.null(manual_remove)){
    graphData <- dplyr::filter(graphData, !(chnm %in% manual_remove))
  }
  
  return(graphData)
}