
# load packages -----------------------------------------------------------

library(terra)
library(tidyverse)
library(here)


# load data ---------------------------------------------------------------

myrcia.files <- list.files(
  here("output", "models_binary_prediction", "raster_binary_05_degree"),
  full.names = T, recursive = T, pattern = ".tif")

myrcia.files <- grep(myrcia.files, pattern = "few_occ|thr_site", value = T)

bin.stk <- rast(myrcia.files)

myrcia_binary_df_05_degree <- bin.stk %>% 
  terra::as.data.frame(xy=T, na.rm = F) %>% 
  pivot_longer(
    cols = 3:309,
    names_to = "species", 
    values_to = "presence") %>% 
  filter(presence == 1) %>% 
  mutate(species = word(species, 1, 2, sep = "_")) %>% 
  arrange(species)


myrcia_binary_df_05_degree %>% 
 write.csv(
   here("output", "models_binary_prediction", "myrcia_binary_df_05_degree.csv"), 
   row.names = F
   )


 