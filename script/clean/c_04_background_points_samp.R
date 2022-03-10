
library(dismo)
library(sp)
library(raster)

load(here("data","occ_myrcia.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))

sp.occ <- 
  occ.myrcia %>% 
  dplyr::select(decimalLongitude, decimalLatitude) %>% 
  SpatialPoints()


bias.std <- readRDS( here("data", "bias_raster_std.rds"))

# ensure that background points have environmental data
env.na <- is.na(values(sum(env.layers)))
values(bias.std)[env.na] <- NA

# sample bg points
bg <- randomPoints(bias.std, n = 20e3, p = sp.occ,  prob = T, excludep = F)

# bias and bg sample
plot(bias.std)
points(bg, pch = 16, cex = .1)

# transform to spatial points and define crs
bg.spatial <- SpatialPoints(bg)
crs(bg.spatial) <- crs(bias.std)

# save
saveRDS(bg.spatial, here("data", "background_pts.rds"))
  



