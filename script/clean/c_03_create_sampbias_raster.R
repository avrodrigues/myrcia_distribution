#require("devtools")
#install_github("azizka/sampbias")
library(sampbias)
library(here)
library(raster)
library(sp)
library(dplyr)

# load env data -----------------------------------------------------------

load(here("data", "occ_df_kew.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))

occ.to.bias <- 
  occ.df.kew %>% 
  dplyr::select(species, decimalLongitude, decimalLatitude)
  

# calculate bias
bias.calc <- 
calculate_bias(x = occ.to.bias, 
               inp_raster = env.layers[[1]])

bias.calc$bias_estimate
# project bias
proj <- project_bias(bias.calc)

# ensemble of bias factors
bias.occ <- proj$cities.roads.rivers

# view bias raster
plot(bias.occ)

# rescale to have the same extent and resoluton of 'env.layers'
bias.resampled <- resample(bias.occ, env.layers, method="bilinear")

# standardize by max value
bias.std <- bias.resampled / max(values(bias.resampled), na.rm = T)
plot(bias.std)


saveRDS(bias.std, here("data", "bias_raster_std.rds"))



