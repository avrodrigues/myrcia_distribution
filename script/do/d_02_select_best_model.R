
# load packages -----------------------------------------------------------

library(ENMeval)
library(tidyverse)
library(raster)
library(here)

# load model results ------------------------------------------------------

l.files.mod <- list.files(
  here("output", "models", "tuned_models"),
  full.names = T
)


# select the best model among tuned models for each species ---------------

## the best model is those which has:
## - lowest AICc (all models with delta AICc < 2 is considered equivalent);
## - highest average of Continuos boyce index;
## - lowest number of coefcients.

# directory to save best model rasters
dir.save <- here("output", "models", "raster_best_models")

# best models evaluation
best.mod.eval <- vector("list", length(l.files.mod))

# models without average of Continuos boyce index (NA)
none.model <- list()
idx <- 1


### select best model, save evaluation and raster
for(i in seq_along(l.files.mod)){
  species <- str_sub(l.files.mod[[i]], 124, nchar(l.files.mod[[i]])-11)
  
  tuned.mod <- readRDS(l.files.mod[i])
  e.res <- eval.results(tuned.mod)
  
  best.mod <- 
    e.res %>% 
    filter(delta.AICc < 2) %>% 
    filter(cbi.val.avg == max(cbi.val.avg, na.rm = T)) %>% 
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

saveRDS(best.mod.eval.df, here("output", "models", "best_models_eval_stats.rds"))
saveRDS(none.model, here("output", "models", "species_none_best_models_eval_stats.rds"))

