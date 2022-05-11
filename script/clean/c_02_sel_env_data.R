# load packages -----------------------------------------------------------

library(sp)
library(raster)
library(here)
library(usdm)
library(tidyverse)

# load env data -----------------------------------------------------------

env.files <- list.files(here("data", "env"), 
                        pattern = "clim|soil",
                        full.names = T)

env.r <- stack(env.files)

# all layers must have values. If any is NA, then all layers must be NA.

#finding NA values and replace
df.env <- as.data.frame(env.r)
not.all.nas <- which(
  rowSums(is.na(df.env)) != 0 &
  rowSums(is.na(df.env[,])) < 9
  )

# replace value by NA when any column in a row have NA
df.env2 <- df.env
for(i in 1:ncol(df.env2)){
  df.env2[not.all.nas, i ] <- NA
}

# set values to raster
env.r2 <- setValues(env.r, as.matrix(df.env2))

# VIF of variables
usdm::vif(env.r2[[-c(2:4)]])

saveRDS(env.r2[[-c(2:4)]], here("data", "env", "env_layers.rds"))
