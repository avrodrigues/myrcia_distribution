
# load packages -----------------------------------------------------------

library(ENMeval)
library(tidyverse)
library(raster)
library(here)

# load model results ------------------------------------------------------

l.files.mod <- list.files(
  here("output", "models", "LOO", "LOO_tuned_models"),
  full.names = T
)[-1]


# select the best model among tuned models for each species ---------------

## the best model is those which has:
## - lowest AICc (all models with delta AICc < 2 is considered equivalent);
## - lowest omision rate 10 percentile
## - highest average of Average AUC;
## - lowest number of coefcients.

# directory to save best model rasters
dir.save <- here("output", "models", "LOO", "raster_best_models")

if(!dir.exists(dir.save)) dir.create(dir.save)

# best models evaluation
best.mod.eval <- vector("list", length(l.files.mod))

# models without average of Continuos boyce index (NA)
none.model <- list()
idx <- 1


### select best model, save evaluation and raster
for(i in seq_along(l.files.mod)){
  species <- str_sub(l.files.mod[[i]], 132, nchar(l.files.mod[[i]])-11)
  
  tuned.mod <- readRDS(l.files.mod[i])
  e.res <- eval.results(tuned.mod)
  
  best.mod <- 
    e.res %>% 
    filter(delta.AICc < 2) %>% 
    filter(or.10p.avg == min(or.10p.avg, na.rm = T)) %>% 
    filter(auc.val.avg == max(auc.val.avg, na.rm = T)) %>% 
    filter(ncoef == min(ncoef)) %>% 
    slice_head() %>% 
    mutate(species = species)
  
  # species' models without any 'cbi.val.avg' value
  if(nrow(best.mod) == 0){
    none.model[[idx]] <-    
      e.res %>% 
      filter(delta.AICc < 2) %>% 
      mutate(species = species)
    
    idx <- idx + 1
    next()
  }
  
  # extract and save raster prediction
  r.pred <- eval.predictions(tuned.mod)[[best.mod$tune.args]]
  
  str.b.mod <- str_replace_all(as.character(best.mod$tune.args), "\\.", "_")
  path.save <- paste0(
    dir.save, "/", species, "_", str.b.mod
  )
  writeRaster(r.pred, path.save, format = "GTiff", overwrite = T)
  
  # save best model evaluation stats
  best.mod.eval[[i]] <- best.mod
  
}

best.mod.eval.df <- bind_rows(best.mod.eval)

saveRDS(best.mod.eval.df, here("output", "models", "LOO", "best_models_eval_stats.rds"))

if(length(none.model) != 0){
  saveRDS(none.model, here("output", "models", "species_none_best_models_eval_stats.rds"))
  
}

