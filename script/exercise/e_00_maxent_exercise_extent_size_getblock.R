###########################################################################
###########################################################################
#
# the aim here is to define a unique size of the buffer used to sample the 
# background points for all species in the study
#
# the main reference for the exercise is the 
# VanDerWal et al. (2009) Ecological Modelling. 


# load packages -----------------------------------------------------------

library(raster)
library(sp)
library(here)
library(dismo)
library(rgeos)
library(sf)
library(ENMeval)
library(dplyr)
library(ROCR)
library(ecospat)
library(ggplot2)


# load data ---------------------------------------------------------------

load(here("data","occ_myrcia.RData"))
env.layers <- readRDS(here("data", "env", "env_layers.rds"))
env.layers.std <- scale(env.layers)
bg.spatial <- readRDS(here("data", "background_pts.rds"))

# sample of species for the exercise --------------------------------------

# species with at least 25 occurrences
spp.names.25 <- 
  occ.myrcia %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  filter(n >= 25) %>% 
  select(species) %>% 
  unlist() %>% 
  as.character()

## species in each group of occ number

spp.strata <- 
  occ.myrcia %>% 
  filter(species %in% spp.names.25) %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  mutate(group = 
           case_when(n >= 25 & n < 50 ~ "1- 25-50",
                     n >= 50 & n < 100 ~ "2- 50-100",
                     n >= 100 & n < 300 ~ "3- 100-300",
                     n >= 300 ~ "4- 300+")
  ) %>% 
  group_by(group) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n), 
         samp.n = round(prop*30, 0))


# sample 30 species to exercise
n.occ.species <- 
  occ.myrcia %>% 
  filter(species %in% spp.names.25) %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  mutate(group = 
           case_when(n >= 25 & n < 50 ~ "1- 25-50",
                     n >= 50 & n < 100 ~ "2- 50-100",
                     n >= 100 & n < 300 ~ "3- 100-300",
                     n >= 300 ~ "4- 300+"),
  )

set.seed(181221)
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

spp.to.exerc <- bind_rows(spp.to.exerc)$species
occ.df.exerc <- occ.myrcia %>% filter(species %in% spp.to.exerc)


# modeling ----------------------------------------------------------------

wgs84 <- crs(env.layers)

eval.list <- vector("list", length(spp.to.exerc))

## loop for modeling
for(i in seq_along(spp.to.exerc)){
  ## select species
  occur <- occ.df.exerc %>% 
    filter(species %in% spp.to.exerc[i]) %>% 
    select(decimalLongitude, decimalLatitude)
  
  occur.env <- raster::extract(env.layers.std, occur, df = TRUE)
  
  
  ## Buffer sizes
  
  buffer.sizes <- seq(100e3, 1000e3, by = 100e3)
  
  sp.occ <- SpatialPoints(occur, proj4string = wgs84)
  
  l.buffer.area <- lapply(buffer.sizes, function(size){
    raster::buffer(sp.occ, size)
  })
  
  eval.df <- data.frame(
    auc = rep(NA, length(buffer.sizes)),
    boyce = rep(NA, length(buffer.sizes)), 
    npres = rep(NA, length(buffer.sizes)),
    buff.sz = rep(NA, length(buffer.sizes)),
    species = rep(NA, length(buffer.sizes))
  )
  
  
  for(b in seq_along(l.buffer.area)){
    
    clima.buffer <- mask(env.layers.std, l.buffer.area[[b]])
    cell.number <- rownames(na.omit(as.data.frame(values(clima.buffer))))
    back.buff <- gIntersection(bg.spatial, l.buffer.area[[b]]) 
    
    ## Pontos de presença e background
    pres <- as.data.frame(sp.occ) 
    names(pres) <- c("x", "y")
    back <- as.data.frame(back.buff)
    
    Pback <- rbind(pres, back)
    Pback$Species <- c(rep(1, nrow(pres)),
                       rep(0, nrow(back)))
    
    pa_data <- 
      st_as_sf(Pback, coords = c("x", "y"), 
               crs = wgs84 )
    
    
    ## spatial block
    # spatial blocking by specified range with random assignment
    block <- get.block(pres, back, orientation = "lat_lon")
    
    # extract the raster values for the species points as a dataframe
    mydata <- raster::extract(clima.buffer, pa_data) %>% 
      as.data.frame()
    # adding species column to the dataframe
    mydata$Species <- pa_data$Species
    
    
    # extract the foldIDs in SpatialBlock object 
    # created in the previous section
    # the folds (list) works for all three blocking strategies
    mydata$folds <- do.call(c, block)
    
    # create a data.frame to store the prediction of each fold (record)
    testTable <- mydata %>% dplyr::select(Species, folds)
    testTable$pred.raw <- NA
    testTable$pred.cloglog <- NA
    
    for(k in seq_len(max(mydata$folds))){
      
      # extracting the training and testing indices
      # this way works with folds list (but not foldID)
      
      
      trainSet <-  mydata %>% 
        filter(folds != k) # training set indices
      testSet <-  mydata %>% 
        filter(folds == k) # testing set indices
      
      m <- maxent(trainSet[, 1:6], trainSet[, 7])
     
      
      test.rows <- testTable$folds == k
      testTable$pred.raw[test.rows] <- predict(m, 
                                         testSet, 
                                         args=c("outputformat=raw","outputformat=cloglog")) 
      
      testTable$pred.cloglog[test.rows] <- predict(m, 
                                                   testSet, 
                                                   args=c("outputformat=cloglog")) 
     
      
      
      
    }
    
    
    pred.auc <- prediction(testTable$pred.raw, testTable$Species)
    auc <- performance(pred.auc,'auc')@y.values[[1]]
    
    
    pred.boyce <- 
      testTable %>% 
      filter(Species == 1) %>% 
      dplyr::select(pred.cloglog)
    
    boyce <- ecospat.boyce(testTable$pred.cloglog, pred.boyce, PEplot = T)$Spearman.cor
    
    #maxent.model[[k]] <- m
    
    eval.df$auc[b] <- auc
    eval.df$boyce[b] <- boyce
    eval.df$npres[b] <- nrow(pred.boyce)
    eval.df$buff.sz[b] <- buffer.sizes[b]/1000
    eval.df$species[b] <- spp.to.exerc[i]
    
    
    
  }
  eval.list[[i]] <- eval.df
  cat(paste(i, "of", length(spp.to.exerc)), fill = T)
}


# plot summarized results -------------------------------------------------
library(patchwork)


eval.df <- bind_rows(eval.list) 

summ.result <- 
  eval.df %>% 
  group_by(buff.sz) %>% 
  summarise(auc.mean = mean(auc, na.rm = T),
            auc.se = sqrt(sd(auc, na.rm = T)),
            boyce.mean = mean(boyce, na.rm = T),
            boyce.se = sqrt(sd(boyce, na.rm = T))) 

bind_rows(eval.list) %>% 
ggplot(aes(x = buff.sz, y = boyce, group = species)) +
  geom_line(size = 0.2) +
  ylim(-1,1)



leg_lines <- 
  tribble(
    ~x, ~y, ~type, ~size,
    105, -.7, "species", 1,
    170, -.7, "species", 1,
    105, -.8, "average", .3,
    170, -.8, "average", .3
  ) 

leg_label <- 
  tribble(
    ~x, ~y, ~label, 
    180, -.7, "Species", 
    180, -.8, "Average", 
  )


(eval.res.plot <- 
summ.result %>% 
    ggplot() +
    geom_line(
      data = eval.df, 
      aes(x = buff.sz, y = boyce, group = species),
      size = 0.1, linetype = 2) +
    geom_line(aes(x = buff.sz, y = boyce.mean)) +
    geom_point(aes(x = buff.sz, y = boyce.mean)) +
    scale_x_continuous(breaks = seq(100, 1000, by = 100)) +
    labs(y = "continous Boyce Index", 
         x = "Buffer size (km)") +
    ylim(-1,1) +
    geom_line(
      data = leg_lines,
      aes(x = x, y = y, linetype = type), 
      size = leg_lines$size,
      ) +
    geom_text(
      data = leg_label, 
      aes(x = x, y = y, label = label),
      hjust = 0
    ) +
   theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.2),
      legend.position = 'none'
      )
  )


ggsave(
  here::here("output/fig/00_buffer_size_eval_ENMeval.png"),
  eval.res.plot,
  height = 5,
  width = 7,
  units = "in", 
  device = "png"
)



