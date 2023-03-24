## downscale species binary prediction to 0.5 degree
library(raster)
library(here)
library(stringr)


myrcia.dir <- list.files(
  here("output", "models", "raster_binary"),
  full.names = T)

save_dirs <- here(
  "output", "models", "raster_binary_05_degree",
  c("few_occ", "thr_site")
  )

for(path in save_dirs) {
  if(!dir.exists(path)) dir.create(path)
}

for(k in seq_along(myrcia.dir)){
  
  file_from <- list.files(
    here(myrcia.dir[k]),
    full.names = T
  )
  
  dir_to <- save_dirs[k]
  
  for(i in seq_along(file_from)){
    
    bin.05 <- raster::aggregate(
      raster::raster(file_from[i]), 
      3,
      fun = max
      )
    
    file.name <- str_split(file_from[i], "./")[[1]]
    file.name <- str_remove(file.name[length(file.name)], ".tif")
    
    path.save <- paste0(
      dir_to, "/", file.name
    )
    
    writeRaster(bin.05, path.save, format = "GTiff", overwrite = T)
  }

}



