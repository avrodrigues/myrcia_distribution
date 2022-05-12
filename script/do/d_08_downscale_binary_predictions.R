## downscale species binary prediction to 0.5 degree
library(raster)
library(here)
library(stringr)


myrcia.bin <- list.files(here("output", "models", "raster_binary_prediction"),
                         full.names = T)

bin_05_degree <- here("output", "models", "raster_binary_prediction", "bin_05_degree")
if(!dir.exists(bin_05_degree)) dir.create(bin_05_degree)



for(i in seq_along(myrcia.bin)){
  bin.05 <- raster::aggregate(raster(myrcia.bin[i]), 3, fun = max)
  file.name <- str_split(myrcia.bin[i], "./")[[1]]
  file.name <- str_remove(file.name[length(file.name)], ".tif")
  
  path.save <- paste0(
    bin_05_degree, "/", file.name
  )
  writeRaster(bin.05, path.save, format = "GTiff", overwrite = T)
  
}
myrcia.bin[i]

