#### Create binary predictions for each species
#### Method: pS-SDM + PRR - Scherrer et al 2018. Methods in Ecol. Evol. 
#### pS-SDM + PRR is a Site-level threshold
#### sum the probability of presence of all species in each site (pS-SDM)
#### then pick the species with the highest prob. until the richeness is equal 
#### to the sum of pS-SDM. 


# load packages -----------------------------------------------------------

library(raster)
library(here)
library(tidyverse)


# Load data ---------------------------------------------------------------

myrcia.good <- list.files(
  here("output", "models", "CV_spatial_block", "raster_good_models"),
  full.names = T
  )

myrcia.good.LOO <- list.files(
  here("output", "models", "CV_leave_one_out", "raster_good_models"), 
  full.names = T
  )


# stacked SDM -------------------------------------------------------------

myrcia.stk <- raster::stack(c(myrcia.good, myrcia.good.LOO))
myrcia.stk.df <- as.data.frame(myrcia.stk, xy = T)

# richness by sites
myrcia.rich.vec <- myrcia.stk.df %>% 
  dplyr::select(-c(1:2)) %>% 
  mutate(rich = rowSums(., na.rm = T)) %>% 
  pull(rich)

names(myrcia.stk.df)

myrcia.df.long <- 
  myrcia.stk.df %>% 
    mutate(site = 1:nrow(myrcia.stk.df), .before = x) %>% 
  pivot_longer(
    4:last_col(), 
    names_to = "species", 
    values_to = "prob"
  ) %>% 
  drop_na(prob)


myrcia.site.sel <- tibble(
  site = integer(),
  x = numeric(),
  y = numeric(),
  species = character(),
  prob = numeric(),
  presence = integer()
)

site.vec <- which(round(myrcia.rich.vec) > 0)

for(i in site.vec){

    add.to.df <- 
      myrcia.df.long %>% 
      filter(site == i) %>% 
      arrange(desc(prob)) %>% 
      slice_head(n = round(myrcia.rich.vec[i])) %>% 
      mutate(presence = 1)
    
    myrcia.site.sel <- 
      myrcia.site.sel %>% 
      add_row(add.to.df)
  
  
  cat(paste("site", i), fill = T)
} # end in 226653



# save binary predictions -------------------------------------------------
dir.save <- here(
  "output", "models_binary_prediction", "raster_binary_10_minutes", "thr_site"
  )

if(!dir.exists(dir.save)) dir.create(dir.save)

sp.file.names <- unique(myrcia.site.sel$species)
xy.site <- xyFromCell(myrcia.stk, myrcia.site.sel$site) 


for(i in seq_along(sp.file.names)){
  sp.xyz <- 
    bind_cols(myrcia.site.sel, xy.site) %>% 
    filter(species == sp.file.names[i])
  
  r.generic <- myrcia.stk[[sp.file.names[i]]]
  values.vec <- integer(length = ncell(r.generic))
  values.vec[sp.xyz$site] <- 1
  r.generic <- raster::setValues(r.generic, values.vec) 
  
  path.save <- paste0(
    dir.save, "/", sp.file.names[i]
  )
  writeRaster(r.generic, path.save, format = "GTiff", overwrite = T)
}




