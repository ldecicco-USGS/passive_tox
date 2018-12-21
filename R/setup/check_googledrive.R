# Modified from:
# https://github.com/ropensci/drake/issues/252

drive_get_datetime_modified <- function(dr_id){

  datetime_modified <- dr_id %>% 
    drive_get %>% 
    mutate(modified = lubridate::as_datetime(purrr::map_chr(drive_resource, "modifiedTime"))) %>% 
    pull
  
  return(datetime_modified)
}

drive_download_gd <- function(dr_id, path, time_stamp){
  
  drive_download(dr_id, path, overwrite = TRUE)
  
  if(file.exists(path)){
    return(path)
  } else {
    return(NULL)
  }
  
}
