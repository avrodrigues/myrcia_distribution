# load packages -----------------------------------------------------------

library(ENMeval)
library(tidyverse)
library(raster)
library(here)


# which species to re-do ------------------------------------------------------

best.mod.eval.df <- readRDS(here("output", "models", "best_models_eval_stats.rds"))


species.worst.mods <- 
  best.mod.eval.df %>% 
  as_tibble() %>% 
  filter(cbi.val.avg < 0.5) 

redo.species <- 
  species.worst.mods %>% 
  select(species) %>% 
  unlist() %>% 
  setNames(NULL) %>% 
  str_replace("_", " ") %>% 
  str_to_sentence()
  

# load data ---------------------------------------------------------------

load(here("data","occ_myrcia.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))
wgs84 <- crs(env.layers)
bg.spatial <- readRDS(here("data", "background_pts.rds"))


# modeling ----------------------------------------------------------------
library(raster)
library(sp)
library(here)
library(dismo)
library(rgeos)
library(sf)
library(ENMeval)
library(kuenm)
library(dplyr)
library(ROCR)
library(ecospat)


#dir.save <- here("output", "models", "tuned_models")
tuned.mod <- list(length(redo.species))
for(i in seq_along(redo.species)){
  
  occ <- 
    occ.myrcia %>% 
    filter(species %in% redo.species[i]) 
  
  env.data <- raster::extract(
    env.layers,
    occ[,c("decimalLongitude", "decimalLatitude")]
  ) %>% as.data.frame()
  
  grid.res <-
    sapply(env.data, range, na.rm = T) %>% 
    t() %>% 
    apply(1, diff)/8 
  
  grid.res <- round(grid.res)
  
  cat(
    paste(
      "conducting env filtering - species", i, "of",  length(redo.species)
      ),
    fill = T
  )
  occur <- 
    env_grid_filter(
      occ, 
      env.data,
      grid.res
    ) %>% 
    select(decimalLongitude, decimalLatitude)
  
  ## Buffer extent
  sp.occ <- SpatialPoints(occur, proj4string = wgs84)
  buffer.area <-  raster::buffer(sp.occ, 500e3)
  
  clima.buffer <- raster::mask(env.layers, buffer.area)
  #cell.number <- rownames(na.omit(as.data.frame(values(clima.buffer))))
  back.buff <- gIntersection(bg.spatial, buffer.area) 
  
  ## Pontos de presença e background
  pres <- as.data.frame(sp.occ) 
  names(pres) <- c("x", "y")
  back <- as.data.frame(back.buff)
  
  cat(
    paste(
      "tuning SDM - species", i, "of",  length(redo.species)
    ),
    fill = T
  )
  
  tuned.mod[[i]] <- ENMevaluate(occs = pres, envs = clima.buffer, bg = back, 
                           algorithm = 'maxent.jar', partitions = 'block',
                           tune.args = list(
                             fc = c('L','Q', 'LQ','LQP','LQH','LQHP'), 
                             rm = seq(0.5, 5, by = 0.5)), 
                           parallel = T, numCores = 3, 
                           quiet = T)
  
 #species <-  spp.names.25[i] %>% 
 #  stringr::str_replace(" ", "_") %>% 
 #  stringr::str_to_lower()
 #
 #path <- paste0(dir.save, "/", species, "_maxent.rds")
 #
 #saveRDS(tuned.mod, path)
  
  cat(paste("species", i, "of", length(redo.species), "modeled"), fill = T)
}



# compare model fit  ------------------------------------------------------

#modeling without env fintering
best.mod.eval.df <- readRDS(
  here("output", "models", "best_models_eval_stats.rds")
)

# model with env filtering
redo.species_ <- str_replace_all(redo.species, " ", "_") %>% 
  str_to_lower()

length(tuned.mod)
names(tuned.mod) <- redo.species_

best.mod.redone <- list(length(tuned.mod))

for(i in seq_along(tuned.mod)){
  best.mod.redone[[i]]<- 
  eval.results(tuned.mod[[i]]) %>% 
    filter(delta.AICc < 2) %>% 
    filter(cbi.val.avg == max(cbi.val.avg, na.rm = T)) %>% 
    filter(or.10p.avg == min(or.10p.avg, na.rm = T)) %>% 
    filter(ncoef == min(ncoef)) %>% 
    slice_head() %>% 
    mutate(species = redo.species_[i])
}


col.select <- c(
  "species",
  "tune.args", 
  "auc.val.avg",
  "cbi.val.avg", 
  "or.10p.avg", 
  "w.AIC",
  "ncoef"
)

# comparing
best.mod.redone.eval.df <- 
  bind_rows(best.mod.redone) %>% 
    select(all_of(col.select)) %>% 
    mutate(env_filter = TRUE)

best.mod.eval.df.original <- 
  best.mod.eval.df %>% 
  filter(species %in% best.mod.redone.eval.df$species) %>% 
  select(all_of(col.select)) %>% 
  mutate(env_filter = FALSE)

compare.best.mod.df <- bind_rows(
  best.mod.redone.eval.df, best.mod.eval.df.original
) %>% 
  arrange(species)
blue_gold <- c(
  "#3B3C68",
  "#C1DBCB",
  "#DFA42F",
  "#846435"
)

blue_gold_red <- c(
  "#303260",
  "#88beca",
  "#e4b434",
  "#7b5a28",
  "#e19296",
  "#c62f22"
  
)

ggplot(compare.best.mod.df, aes(y = as.factor(species), x = or.10p.avg)) +
  geom_line(aes(group = species), size = 2, color = greys_10[6]) +
  geom_point(aes(color = env_filter), size = 4) +
  scale_color_manual(values = blue_gold_red[c(1,2)]) +
  theme_bw()

