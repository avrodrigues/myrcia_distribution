
# load packages -----------------------------------------------------------

library(terra)
library(tidyverse)
library(here)


# load data ---------------------------------------------------------------

bin_05_degree <- here("output", "models", "raster_binary_prediction", "bin_05_degree")
files_bin_05_degree <- list.files(bin_05_degree, full.names = T)


bin.stk <- rast(files_bin_05_degree)

myrcia_binary_df_05_degree <- 
bin.stk %>% 
  as.data.frame(xy=T) %>% 
  pivot_longer(
    cols = myrcia_aethusa_fc_LQP_rm_1:myrcia_zuzygium_fc_LQP_rm_2,
    names_to = "species", 
    values_to = "presence") %>% 
  filter(presence == 1) %>% 
  mutate(species = word(species, 1, 2, sep = "_")) %>% 
  arrange(species)



myrcia_binary_df_05_degree %>% 
 saveRDS(here("output", "models", "myrcia_binary_df_05_degree.rds"))


 