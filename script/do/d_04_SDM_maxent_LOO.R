library(raster)
library(sp)
library(here)
library(dismo)
library(rgeos)
library(sf)
library(ENMeval)
library(kuenm)
library(tidyverse)
library(ROCR)
library(ecospat)

# load data ---------------------------------------------------------------

load(here("data","occ_myrcia.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))
wgs84 <- crs(env.layers)
bg.spatial <- readRDS(here("data", "background_pts.rds"))

# species names -----------------------------------------------------------

# species with at number of occ between 6  and 14
spp.names.6.24 <- 
  occ.myrcia %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  filter(n > 5 & n < 25) %>% 
  dplyr::select(species) %>% 
  unlist() %>% 
  as.character()



# modeling ----------------------------------------------------------------
dir.save <- here("output", "models", "CV_leave_one_out", "LOO_tunned_models")

for(i in seq_along(spp.names.6.24)){
  ## select species
  occur <- occ.myrcia %>% 
    filter(species %in% spp.names.6.24[i]) %>% 
    select(decimalLongitude, decimalLatitude)
  
  ## Buffer extent
  sp.occ <- SpatialPoints(occur, proj4string = wgs84)
  buffer.area <-  raster::buffer(sp.occ, 700e3)
  
  clima.buffer <- mask(env.layers, buffer.area)
  #cell.number <- rownames(na.omit(as.data.frame(values(clima.buffer))))
  back.buff <- gIntersection(bg.spatial, buffer.area) 
  
  ## Pontos de presença e background
  pres <- as.data.frame(sp.occ) 
  names(pres) <- c("x", "y")
  back <- as.data.frame(back.buff)
  
  tuned.mod <- ENMevaluate(occs = pres, envs = clima.buffer, bg = back, 
                           algorithm = 'maxent.jar', partitions = 'jackknife',
                           tune.args = list(
                             fc = c('L','Q', 'LQ','LQP','LQH','LQHP'), 
                             rm = seq(0.5, 5, by = 0.5)), 
                           parallel = T, numCores = 3, 
                           quiet = T)
  
  species <-  spp.names.6.24[i] %>% 
    stringr::str_replace(" ", "_") %>% 
    stringr::str_to_lower()
  
  path <- paste0(dir.save, "/", species, "_maxent.rds")
  
  saveRDS(tuned.mod, path)
}

