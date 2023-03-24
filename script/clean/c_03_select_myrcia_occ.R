library(tidyverse)
library(naturaList)
library(here)
source(here("function", "env_grid_filter.R"))

# taxonomic revision made based on POWO and filter for occ recorded after 1979
load(here::here("data", "occ_df_kew.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))

# select Myrcia occurrences
occ.myrcia.raw <- 
occ.df.kew %>% 
  filter(str_detect(species, "Myrcia "))

n_sp_raw <- occ.myrcia.raw$species %>% unique() %>% length()

## Filter occurrence records in environmental space with the best classification 
## of species ID

## #######################
## create grid bins for all region!!!
## #########################################

env.data.neo <- raster::extract(
  env.layers,
  occ.myrcia.raw[,c("decimalLongitude", "decimalLatitude")]
) %>% as.data.frame()

grid.res <-
  sapply(env.data.neo, range, na.rm = T) %>% 
  t() %>% 
  apply(1, diff)/10 

grid.res <- round(grid.res)



sp.names <- unique(occ.myrcia.raw$species)
g.filter <- vector("list", length(sp.names))

for(i in seq_along(sp.names)){
  occ.temp <- occ.myrcia.raw %>% 
    filter(species %in% sp.names[i])

  env.data <- raster::extract(
    env.layers,
    occ.temp[,c("decimalLongitude", "decimalLatitude")]
  ) %>% as.data.frame()
  

  
  if(nrow(occ.temp) > 1){
    
    g.filter[[i]] <-   env_grid_filter(
      occ.temp, 
      env.data,
      grid.res
    ) 
  }else{ 
    g.filter[[i]] <- occ.temp
    next() 
    }
  
  
}
  
occ.myrcia <- 
  bind_rows(g.filter)

nrow(occ.myrcia)
nrow(occ.myrcia.raw)

plot(env.layers$clim_var_bio_4)
points( occ.myrcia.raw[,c("decimalLongitude", "decimalLatitude")], cex = 0.2)


plot(env.layers$clim_var_bio_4)
points( occ.myrcia[,c("decimalLongitude", "decimalLatitude")], cex = 0.2)


# save
save(occ.myrcia, file = here::here("data", "occ_myrcia.RData"))

