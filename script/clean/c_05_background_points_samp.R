

# load packages -----------------------------------------------------------


library(dismo)
library(sp)
library(raster)


# load data ---------------------------------------------------------------


env.layers <- readRDS(here("data", "env", "env_layers.rds"))
wgs84 <- crs(env.layers)
bias.std <- readRDS( here("data", "bias_raster_std.rds"))

# ensure that background points have environmental data
env.na <- is.na(values(sum(env.layers)))
values(bias.std)[env.na] <- NA

# sample bg points
bg <- randomPoints(bias.std, n = 20e3, prob = T, excludep = F)

# bias and bg sample
plot(bias.std)
points(bg, pch = 16, cex = .1)

# transform to spatial points and define crs
bg.spatial <- SpatialPoints(bg)
crs(bg.spatial) <- wgs84

# save
saveRDS(bg.spatial, here("data", "background_pts.rds"))
  



