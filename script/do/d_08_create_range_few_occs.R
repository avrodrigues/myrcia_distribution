# load packages -----------------------------------------------------------

library(raster)
library(sf)
library(here)
library(tidyverse)
library(rnaturalearth)

# load data ---------------------------------------------------------------

load(here("data","occ_myrcia.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))
wgs84 <- crs(env.layers)

# defining geographic region
neo_coastline <- read_sf(
  here("data", "shapes", "Continentes", "level1.shp")
  )

ext <- extent(env.layers)
neo_coastline <- st_crop(neo_coastline, ext) %>% st_union()
sf::st_crs(neo_coastline) <- 4326

# species names -----------------------------------------------------------

# species with at less than 5 occurrences
# and only identified by specialists
spp.names.few.occ <- 
  occ.myrcia %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  filter(n <= 5) %>% 
  dplyr::select(species) %>% 
  unlist() %>% 
  as.character()

dir.save <- here("output", "models_binary_prediction", "raster_binary_10_minutes", "few_occ")

for(i in seq_along(spp.names.few.occ)){
  occ.temp.sf <- 
    occ.myrcia %>% 
    filter(species == spp.names.few.occ[i]) %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude")) 
  
  pt <- st_union(occ.temp.sf)
  
 
  sp.pol <- sf::st_buffer(pt, 0.5)
  
  geo <- sf::st_geometry(sp.pol)
  sf::st_crs(geo) <- 4326
  sp.range <- suppressMessages(sf::st_intersection(geo, neo_coastline))
  
  plot(neo_coastline, main = spp.names.few.occ[i])
  plot(sp.range, col= "#BD0B28", add = T)

  r.range <- 
    terra::vect(sp.range) %>% 
    terra::rasterize(terra::rast(env.layers)) 
  
  sp_name <- str_to_lower(spp.names.few.occ[i]) %>% 
    str_replace_all(" ", "_") %>% 
    paste0("_few_occ.tif")
  
  path.save <- here(dir.save, sp_name)

  
  terra::writeRaster(
    r.range, 
    filename = path.save
  )
}


