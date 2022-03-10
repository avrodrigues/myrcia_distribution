
# load packages -----------------------------------------------------------


library(sf)
library(here)
library(raster)
library(tidyverse)
library(stars)

# load data -----------------------------------------------------------

load(here("data","occ_myrcia.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))
env.stars <- st_as_stars(env.layers$clim_var_bio_1)
world <- read_sf("data/shapes/Continentes/level1.shp")
st_crs(world) <- "WGS84"

# distribution of species with less than 5 occurrences --------------------

# species with < 10 occurrences
spp.less.5 <- 
  occ.myrcia %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  filter(n < 5) %>% 
  dplyr::select(species) %>% 
  unlist() %>% 
  as.character()

# select occurrences
occ.less.5 <- 
  occ.myrcia %>% 
  filter(species %in% spp.less.5)

# transform to sf
occ.less.5.sf <- st_as_sf(
  occ.less.5,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs =  "WGS84"
  )

# create a buffer of 0.5 decimal degree around occurences
occ.less.5.buffer <- st_buffer(occ.less.5.sf, dist = 0.5) 

# directory to save best model rasters
dir.save <- here("output", "models", "raster_few_occs")

for(i in seq_along(spp.less.5)){
  
  # union of buffers
  union.buf <-
  occ.less.5.buffer %>% 
    filter(species == spp.less.5[i]) %>% 
    st_union() 
  
  # create dist raster
  r.dist <- fasterize::fasterize(sf::st_sf(a = 1, union.buf),
                         env.layers[[1]])
    
  # save raster
  species.save <- spp.less.5[i] %>% 
    str_replace_all(" ", "_") %>% 
    tolower()
  
  path.save <- paste0(
    dir.save, "/", species.save, "_buffer"
  )
  writeRaster(r.dist, path.save, format = "GTiff", overwrite = T)
    
}


# distribution of species with 5 to 9 occurrences -------------------------

spp.btwn.5.10 <- 
  occ.myrcia %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  filter(n >= 5 & n < 10) %>% 
  select(species) %>% 
  unlist() %>% 
  as.character()

# select occurrences
occ.btwn.5.10  <- 
  occ.myrcia %>% 
  filter(species %in% spp.btwn.5.10)

# transform to sf
occ.btwn.5.10.sf <- st_as_sf(
  occ.btwn.5.10,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs =  "WGS84"
)

for(i in seq_along(spp.btwn.5.10)){
  # create convex hull and intersects with continental land
  sp.dist.sf <- 
    occ.btwn.5.10.sf %>% 
    filter(species == spp.btwn.5.10[i]) %>% 
    st_union() %>% 
    st_convex_hull() %>% 
    st_buffer(dist = 0.5) %>% 
    st_intersection(world) %>% 
    st_union()
  
  # create dist raster
  r.dist.btwn.5.10 <- fasterize::fasterize(
    sf::st_sf(a = 1, sp.dist.sf),
    env.layers[[1]]) 
  
  # save raster
  species.save <- spp.btwn.5.10[i] %>% 
    str_replace_all(" ", "_") %>% 
    tolower()
  
  path.save <- paste0(
    dir.save, "/", species.save, "_convex_hull"
  )
  
  writeRaster(r.dist.btwn.5.10, path.save, format = "GTiff", overwrite = T)
  
}

