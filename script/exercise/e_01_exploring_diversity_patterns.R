
# load packages -----------------------------------------------------------

library(sf)
#library(stars)
library(tidyverse)
library(here)
library(terra)
library(distances)


# load data ---------------------------------------------------------------


list.good.raster <- list.files(
  here("output", "models", "raster_good_models"), 
  full.names = T
)

list.buff.raster <- list.files(
  here("output", "models", "raster_few_occs"),
  pattern = "buffer", 
  full.names = T
)

list.chul.raster <- list.files(
  here("output", "models", "raster_few_occs"),
  pattern = "convex_hull", 
  full.names = T
)


r.good <- rast(list.good.raster)
r.buff <- rast(list.buff.raster)
r.chul <- rast(list.chul.raster)

r.all <- rast(c(list.good.raster,
                list.buff.raster,
                list.chul.raster))

rich.good <- sum(r.good, na.rm = T)
plot(rich.good)

rich.buff <- sum(r.buff, na.rm = T)
plot(rich.buff)

rich.chul <- sum(r.chul, na.rm = T)
plot(rich.chul)

rich.total <- sum(r.all, na.rm = T)
plot(rich.total)


##### composition
library(PCPS)
library(vegan)

m.comp <- values(r.all)

rich.df <- 
values(rich.total) %>% 
  as.data.frame() 

cell.na <- is.na(rich.df$sum)

comp <- m.comp[!cell.na, ]
comp[is.na(comp)] <- 0


less.than.3 <- rowSums(comp) < 3.5

less.than.2sp <- colSums(comp[!less.than.3, ]) <= 2.5

dim(comp[!less.than.3, !less.than.2sp])
dim(comp)
comp.hell <- decostand(comp[!less.than.3, !less.than.2sp], "hel")
comp.dist <- distances(comp.hell)

comp.dist <- dist(comp.hell)
dim(comp.dist)

pcoa.myrcia <- ape::pcoa(comp.dist)
