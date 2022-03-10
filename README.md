The goal is to model the species distribution for all *Myrcia* species
with occurrence records available and with good quality.

The modeling strategy is to split the species in two groups, with less
and more than 25 occurrences. The species with few occurrences records
will be modeled using leave-one-out procedure. The other group will be
modeled using spatial blocks.

All models will be fitted with tuned maxent.

## Repository structure

the repository has three main subfolders `script`, `data`, and `output`:

-   **data**
    -   env
-   **output**
    -   models
        -   raster\_best\_models
        -   raster\_good\_models
        -   tuned\_models
    -   fig
-   **script**
    -   clean
    -   do
    -   exercise

### `script` folder

In the `script` folder, the files are named with the following structure
**Number\_LowerCaseLetter\_ScriptName**. - the number indicates the
order that the scripts must be run; - the lower case letter indicates
the type of script and its subfolder (c - clean; d - do; e - exercise) -
the script name has a short description of the task executed in the
script

In the folder `clean` are the scripts to organize the data for the
modeling. In the folder `do` are the scripts for model fitting,
evaluation and model s

d - do

Data has the input data for analysis conducted output has the results of
model fiting and figures script is organized
