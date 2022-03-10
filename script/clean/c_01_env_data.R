
# load packages -----------------------------------------------------------

library(sp)
library(raster)
library(here)
library(rgdal)

# defining region extent
NEO <- list(x = c(-120,-30),
            y = c(-60,35))


# climatic layers ---------------------------------------------------------


l.files.clim <- list.files(here("data", "env", "A_CHELSA_cur_0ka_V1_2B_r10m"), 
                      full.names = T)

bioclim.files <- grep(".tif$", l.files.clim, value = T)

#load layers
r.bioclim <- stack(bioclim.files[c(1,4,7,9,14)])
# crop by neotropical region extent
neo.bioclim <- crop(r.bioclim, NEO)

writeRaster(neo.bioclim, here("data","env","clim_var.tif"), bylayer = T, suffix = "names")

# Soil layers -------------------------------------------------------------


l.files.soil <- list.files(here("data", "env", "solos_HWSD", "data"), 
                      full.names = T)

sel.file.solos <- grep("T_CEC_CLAY|T_CLAY",
                       l.files.soil, value = T)

soil.hwsd <- stack(sel.file.solos)
neo.soil.hwsd <- crop(soil.hwsd, NEO)

neo.soil.hwsd.resampled <- resample(
  neo.soil.hwsd, neo.bioclim[[1]], method="bilinear"
  )

names(neo.soil.hwsd.resampled) <- c("T_CEC_CLAY", "T_CLAY")

writeRaster(neo.soil.hwsd.resampled, here("data","env","soil_var.tif"),
            bylayer = T, suffix = "names")


