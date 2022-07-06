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
  here("output", "models", "raster_good_models"),
  full.names = T
  )

myrcia.good.LOO <- list.files(
  here("output", "models", "LOO", "raster_good_models"), 
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
dir.save <- here("output", "models", "raster_binary", "thr_site")

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



# binary maps fixed threshold ---------------------------------------------
load(here("data", "occ_myrcia.RData"))


myrcia.df.long.2 <- 
myrcia.df.long %>% 
  mutate(species = species %>% 
           word(1, 2, sep = "_") %>% 
           str_to_sentence() %>% 
           str_replace_all("_", " "),
         )  
  
species_modeled <- unique(myrcia.df.long.2$species) %>% 
  str_sort()

occ.sp.xy <- 
occ.myrcia %>% 
  filter(species %in% species_modeled) %>% 
  dplyr::select(species, decimalLongitude, decimalLatitude) 

  
cell.occ <- cellFromXY(
  myrcia.stk[[1]], 
  occ.sp.xy %>% 
    dplyr::select(decimalLongitude, decimalLatitude) %>% 
    set_names(c("x", "y")) %>% 
    as.matrix() 
  )

myrcia.df.to.bin.or10p <- 
occ.sp.xy %>% 
  mutate(
    site = cell.occ, 
    obs_occ = 1
  ) %>% 
  right_join(
    myrcia.df.long.2,
    by = c("species", "site")
  ) %>% 
  as_tibble() %>% 
  select(
    site, x, y, species, prob, obs_occ
  )

threshold.df <- 
myrcia.df.to.bin.or10p %>% 
  filter(obs_occ == 1) %>% 
  group_by(species) %>% 
  summarise(
    thr_min = min(prob), 
    thr_10p = quantile(prob, 0.1)
  ) %>% 
  mutate(species = species  %>% 
           str_to_lower() %>% 
           str_replace_all(" ", "_"))

dir_thr_min <- here("output", "models", "raster_binary", "thr_min")
dir_thr_10p <- here("output", "models", "raster_binary", "thr_10p")

for(i in seq_along(species_modeled)){
  
  sp <- species_modeled[i] %>% 
    str_to_lower() %>% 
    str_replace_all(" ", "_")
  
  thr_min <- threshold.df %>% 
    filter(species == sp) %>%
    pull(thr_min)
  
  thr_10p <- threshold.df %>% 
    filter(species == sp) %>%
    pull(thr_10p)
  
  sp.stk <- str_detect(names(myrcia.stk), sp) %>% which()
  
  bin_thr_min <- myrcia.stk[[sp.stk]] >= thr_min
  bin_thr_10p <- myrcia.stk[[sp.stk]] >= thr_10p
  
  path.save.min <- paste0(
    dir_thr_min, "/", names(myrcia.stk)[sp.stk]
  )
  path.save.10p <- paste0(
    dir_thr_10p, "/", names(myrcia.stk)[sp.stk]
  )
  
  
  writeRaster(bin_thr_min, path.save.min, format = "GTiff", overwrite = T)
  writeRaster(bin_thr_10p, path.save.10p, format = "GTiff", overwrite = T)
}


plot(bin_thr_min)
plot(bin_thr_10p)
