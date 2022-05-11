# load packages -----------------------------------------------------------
if(!require(kuenm)){
  devtools::install_github("marlonecobos/kuenm")
}


library(raster)
library(sp)
library(here)
library(dismo)
library(rgeos)
library(sf)
library(ENMeval)
library(kuenm)
library(tidyverse)
library(ROCR)
library(ecospat)

# load data ---------------------------------------------------------------

load(here("data","occ_myrcia.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))
wgs84 <- crs(env.layers)
bg.spatial <- readRDS(here("data", "background_pts.rds"))

# species names -----------------------------------------------------------

# species with at number of occ between 6  and 14
spp.names.6.24 <- 
  occ.myrcia %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  filter(n > 5 & n < 25) %>% 
  dplyr::select(species) %>% 
  unlist() %>% 
  as.character()

## species in each group of occ number

spp.strata <- 
  occ.myrcia %>% 
  filter(species %in% spp.names.6.24) %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  mutate(group = 
           case_when(n >= 6 & n < 10 ~ "1- 6-9",
                     n >= 10 & n < 15 ~ "2- 10-15",
                     n >= 15 ~ "3- 15-24")
  ) %>% 
  group_by(group) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n), 
         samp.n = round(prop*30, 0))



# sample 30 species to exercise
n.occ.species <- 
  occ.myrcia %>% 
  filter(species %in% spp.names.6.24) %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  mutate(group = 
           case_when(n >= 6 & n < 10 ~ "1- 6-9",
                     n >= 10 & n < 15 ~ "2- 10-15",
                     n >= 15 ~ "3- 15-24")
  ) 

set.seed(300322)
spp.to.exerc <- list()
idx <- 1
for(gr in spp.strata$group){
  n.samples <- filter(spp.strata, group == gr)$samp.n
  
  spp.to.exerc[[idx]] <- 
    n.occ.species %>% 
    filter(group == gr) %>% 
    slice_sample(n = n.samples)
  
  idx <- idx + 1
}



# modeling exercise -------------------------------------------------------


spp.to.exerc <- bind_rows(spp.to.exerc)$species
occ.df.exerc <- occ.myrcia %>% filter(species %in% spp.to.exerc)
dir.save <- here("output", "models", "LOO", "LOO_tuned_models")

buffer_test <- c(500e3, 700e3)
path_test <- paste0(dir.save, "/buffer_size_test/", c("buffer_500", "buffer_700"))

dir.exists(path_test[1])
for(i in seq_along(spp.to.exerc)){
  ## select species
  occur <- occ.df.exerc %>% 
    filter(species %in% spp.to.exerc[i]) %>% 
    dplyr::select(decimalLongitude, decimalLatitude)
  
  for(b in seq_along(buffer_test)){
    ## Buffer extent
    sp.occ <- SpatialPoints(occur, proj4string = wgs84)
    buffer.area <-  raster::buffer(sp.occ, buffer_test[b])
    
    clima.buffer <- mask(env.layers, buffer.area)
    #cell.number <- rownames(na.omit(as.data.frame(values(clima.buffer))))
    back.buff <- gIntersection(bg.spatial, buffer.area) 
    
    ## Pontos de presença e background
    pres <- as.data.frame(sp.occ) 
    names(pres) <- c("x", "y")
    back <- as.data.frame(back.buff)
    
    tuned.mod <- ENMevaluate(occs = pres, envs = clima.buffer, bg = back, 
                             algorithm = 'maxent.jar', partitions = 'jackknife',
                             tune.args = list(
                               fc = c('L','Q', 'LQ','LQP','LQH','LQHP'), 
                               rm = seq(0.5, 5, by = 0.5)), 
                             parallel = T, numCores = 3, 
                             quiet = T)
    
    species <-  spp.to.exerc[i] %>% 
      stringr::str_replace(" ", "_") %>% 
      stringr::str_to_lower()
    
    path <- paste0(path_test[b], "/", species, "_maxent.rds")
    
    saveRDS(tuned.mod, path)
  }
}


# Evaluate models and difference between buffer sizes ---------------------

# file path for model evaluation
l.files.b300 <- list.files(path_test[1], full.names = T)
l.files.b600 <- list.files(path_test[2], full.names = T)

l.files.mod <- list(l.files.b300, l.files.b600)

# best models evaluation
best.mod.eval <- vector("list", length(l.files.mod))
best.mod.eval$b300 <-  vector("list", length(l.files.b300))
best.mod.eval$b600 <-  vector("list", length(l.files.b600))


# models without average of Continuos boyce index (NA)
none.model <- list()
idx <- 1

for(i in seq_along(l.files.mod)){
  for(b in seq_along(l.files.mod[[i]])){
    species <- str_sub(
      l.files.mod[[i]][b],
      160, 
      nchar(l.files.mod[[i]][b])-11
      )
    
    tuned.mod <- readRDS(l.files.mod[[i]][b])
    e.res <- eval.results(tuned.mod)
    
    best.mod <- 
      e.res %>% 
      filter(delta.AICc < 2) %>% 
      #filter(cbi.val.avg == max(cbi.val.avg, na.rm = T)) %>% 
      filter(auc.val.avg == max(auc.val.avg, na.rm = T)) %>% 
      filter(or.mtp.avg == min(or.mtp.avg, na.rm = T)) %>% 
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
    best.mod.eval[[i]][[b]] <- best.mod
  }

  
}

best.mod.eval.df <- bind_rows(best.mod.eval) %>% 
  mutate(buffer.size = rep(c("500", "700"), each = 30)) 
  filter(ncoef != 0)
best.mod.eval.df 


colors.gr <- c("#be4766", "#8ad0d9")

(plot.or.mtp <- 
ggplot(best.mod.eval.df) +
  geom_boxplot(aes(y = or.mtp.avg, group = buffer.size, fill = buffer.size)) +
  scale_fill_manual(values = colors.gr, name = "Buffer size (Km)") +
  theme_bw() +
  scale_y_continuous(name = "Omission rate (mtp)", limits = c(0, 1)) +
  theme(
    legend.position = "bottom", 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) 
)

(plot.auc <- 
  ggplot(best.mod.eval.df) +
  geom_boxplot(aes(y = auc.val.avg, group = buffer.size, fill = buffer.size)) +
  scale_fill_manual(values = colors.gr, name = "Buffer size (Km)") +
  theme_bw() +
  scale_y_continuous(name = "Average AUC", limits = c(0, 1)) +
  theme(
    legend.position = "bottom", 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) 
)


ggplot(best.mod.eval.df) +
  geom_boxplot(aes(y = auc.val.avg, group = buffer.size, fill = buffer.size)) +
  scale_fill_manual(values = colors.gr) +
  theme_bw()


library(patchwork)

plot.or.mtp + plot.auc

auc.t.test <- t.test(auc.val.avg ~ buffer.size, 
                     data = best.mod.eval.df,
                     paired = T )

auc.t.test

or.mtp.t.test <- t.test(or.mtp.avg ~ buffer.size, 
                     data = best.mod.eval.df,
                     paired = T )

or.mtp.t.test

w.AIC.t.test <- t.test(w.AIC ~ buffer.size, 
                        data = best.mod.eval.df,
                        paired = T )

w.AIC.t.test


#(plot.or.mtp.paired <- 
  ggplot(best.mod.eval.df) +
  geom_line(aes(y = or.mtp.avg, x = buffer.size, group = species)) 
  scale_fill_manual(values = colors.gr, name = "Buffer size (Km)") +
  theme_bw() +
  scale_y_continuous(name = "Omission rate (mtp)", limits = c(0, 1)) +
  theme(
    legend.position = "bottom", 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) 
  
  
  best.mod.eval.df %>% 
    group_by(buffer.size, ncoef) %>% 
    summarise(n = n())
  
  
  
  
  