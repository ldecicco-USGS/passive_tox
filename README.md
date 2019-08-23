# Passive Sampler analysis

We're trying the `drake` R package for orchestrating the workflow. 

To initially get the datasets and clean up the data, run the `passive_data_setup.R` script. If there are any issues, you can try to diagnose the problems using `vis_drake_graph`.

The next drake workplan is in "passive_analysis.R".

Once that file is run sucessfully, there will be an Excel file: "data/clean/passive.xlsx". This file can be used for regular `toxEval` workflows.

Updating the data on the Google Drive must follow this pattern to work correctly:
https://www.labnol.org/internet/update-files-in-google-drive/28928/

## Disclaimer

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.