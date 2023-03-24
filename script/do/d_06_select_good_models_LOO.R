# load packages -----------------------------------------------------------

library(ENMeval)
library(tidyverse)
library(raster)
library(here)


# select good models ------------------------------------------------------

best.mod.eval.df <- readRDS(here("output", "models", "CV_leave_one_out", "best_models_eval_stats.rds"))


species.good.mods <- 
best.mod.eval.df %>% 
  as_tibble() %>% 
  filter(auc.val.avg >= 0.70 & or.10p.avg <= 0.25) %>% 
  dplyr::select(species)


dir.best <-here("output", "models", "CV_leave_one_out", "raster_best_models")
dir.good <- here("output", "models", "CV_leave_one_out", "raster_good_models")

if(!dir.exists(dir.good)) dir.create(dir.good)

## copy raster good models to new directory
for(i in species.good.mods$species){
  from.path <- list.files(dir.best, pattern = i, full.names = T)
  to.path <- str_replace(from.path, "raster_best_models", "raster_good_models")
  
  file.copy(from.path, to.path)
}
  
good.models.eval.stats <- 
  best.mod.eval.df %>% 
  filter(auc.val.avg >= 0.70 & or.10p.avg <= 0.25) 

saveRDS(good.models.eval.stats, here("output", "models", "CV_leave_one_out", "good_models_eval_stats.rds"))
