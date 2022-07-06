
# load packages -----------------------------------------------------------

library(terra)
library(tidyverse)
library(here)


# load data ---------------------------------------------------------------

myrcia.dir <- list.files(
  here("output", "models", "raster_binary_05_degree"),
  full.names = T)


myrcia.dir

l.files.to.stk <- lapply(2:4, function(i){
  dir.a <- list.files(myrcia.dir[1], full.names = T)
  dir.b <- list.files(myrcia.dir[i], full.names = T)
  
  c(dir.a, dir.b)
})
  

l.bin.stk <- lapply(l.files.to.stk, rast)

l_myrcia_binary_df_05_degree <- lapply(
  l.bin.stk, function(bin.stk) {
    bin.stk %>% 
      terra::as.data.frame(xy=T, na.rm = F) %>% 
      pivot_longer(
        cols = 3:309,
        names_to = "species", 
        values_to = "presence") %>% 
      filter(presence == 1) %>% 
      mutate(species = word(species, 1, 2, sep = "_")) %>% 
      arrange(species)
  }
)

names(l_myrcia_binary_df_05_degree) <- c(
  "few_occ_thr_10p", 
  "few_occ_thr_min", 
  "few_occ_thr_site"
)

l_myrcia_binary_df_05_degree %>% 
 saveRDS(here("output", "models", "list_myrcia_binary_df_05_degree.rds"))


 