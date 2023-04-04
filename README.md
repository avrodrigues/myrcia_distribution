# Modeling distribution of *Myrcia* species

The goal is to model the species distribution for all *Myrcia* species
with occurrence records available and with good quality.

The modeling strategy was to split the species in two groups, with less
and more than 25 occurrences. The species with few occurrence records
were modeled using leave-one-out cross-validation The other group were
modeled using spatial blocks cross-validation.

All models were fitted with tuned Maxent.

## Repository structure

This repository is organized in four directories:

`data` has the data sets used in the analysis  
`scripts` has the scripts to make the cleaning, analysis and figures  
`output` has the output from analysis steps and also the figures  
`function` has the functions created for the analysis here

The `output` folder has the models results and needs a more detailed
description. In `output/models` you find the models outputs for each
species in `.rds` files inside the
`output/models/Cv_*/*_tunned_models`.  
In `output/models_binary_prediction` you find the raster files (`.tif`)
for prediction in two resolutions, 10 minutes (the scale of the
modeling) and 0.5 degrees (a coarser resolution used to derive diversity
patterns and analysis). You can also find a `.csv` file with species
name and geographical coordinates based on the coarser resolution
raster.
