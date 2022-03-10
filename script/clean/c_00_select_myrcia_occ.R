library(tidyverse)
library(naturaList)
library(here)


# taxonomic revision made based on POWO and filter for occ recorded after 1979
load(here::here("data", "occ_df_kew.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))

# select Myrcia occurrences
occ.myrcia.raw <- 
occ.df.kew %>% 
  filter(str_detect(species, "Myrcia "))

# keep only unique spatial occurrence with the best classification of species ID
sp.names <- unique(occ.myrcia.raw$species)


g.filter <- vector("list", length(sp.names))

for(i in seq_along(sp.names)){
  occ.temp <- occ.myrcia.raw %>% 
    filter(species %in% sp.names[i])
  
  if(nrow(occ.temp) > 1){
    g.filter[[i]] <- grid_filter(occ.temp, r = env.layers$clim_var_bio_1)
  }else{ 
    g.filter[[i]] <- occ.temp
    next() 
    }
  
  
}
  
occ.myrcia <- 
  bind_rows(g.filter)

# save
save(occ.myrcia, file = here::here("data", "occ_myrcia.RData"))
