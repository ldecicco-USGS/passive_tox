# Passive Sampler analysis

We're trying the `drake` R package for orchestrating the workflow. 

To initially get the datasets and clean up the data, run the `passive_data_setup.R` script. If there are any issues, you can try to diagnose the problems using `vis_drake_graph`.

## Converting to OneDrive

Halfway through this project, we switched from GoogleDrive to OneDrive. We're trying to take advantage of the automatic snycing of the data folders. To do this, each collaborator has to:

1. Go to the canonical shared drive location
2. Click the "Sync" button and open in Microsoft OneDrive
3. In Windows Explore* right-click on the shared folder and choose "Alway Keep on this device". (Note: later when you want to remove this, you'll need to right-click -> Settings -> Stop Syncing before you can delete it)
4. Create a system variable with the path to the local copy of the canonical shared drive:

```
rprofile_path = file.path(Sys.getenv("HOME"), ".Rprofile")

write('Sys.setenv(PASSIVE_PATH = "C:/Users/ldecicco/DOI/Corsi, Steven R - Manuscript")',
      rprofile_path, 
      append =  TRUE)

cat('Your Rprofile has been updated to include PASSIVE_PATH
    Please restart R for changes to take effect.')
```
You MUST restart your R session (restart RStudio)!

Then, you can call the chemicalSummary via:
```
chemicalSummary <- readRDS(file = file.path(Sys.getenv("PASSIVE_PATH"),"data","data_for_git_repo","clean","chemical_summary.rds"))
```

If you've run the drake plans, this will also work:
```
loadd(chemicalSummary)
```


## Old GoogleDrive stuff:

The next drake workplan is in "passive_analysis.R".

Once that file is run sucessfully, there will be an Excel file: "data/clean/passive.xlsx". This file can be used for regular `toxEval` workflows.

Updating the data on the Google Drive must follow this pattern to work correctly:
https://www.labnol.org/internet/update-files-in-google-drive/28928/

## Disclaimer

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.