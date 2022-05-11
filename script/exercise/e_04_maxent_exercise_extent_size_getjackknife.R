###########################################################################
###########################################################################
#
# the aim here is to define a unique size of the buffer used to sample the 
# background points for all species in the study with more than 5 occurrences
# and less than 25 occurrences
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
  
  occur.df <- bind_cols(occur, occur.env) %>% 
    tidyr::drop_na() %>% 
    dplyr::select(decimalLongitude, decimalLatitude) 
  
  ## Buffer sizes
  
  buffer.sizes <- seq(100e3, 1000e3, by = 100e3)
  
  sp.occ <- SpatialPoints(occur.df, proj4string = wgs84)
  
  l.buffer.area <- lapply(buffer.sizes, function(size){
    raster::buffer(sp.occ, size)
  })
  
  eval.df <- data.frame(
    avg.auc = rep(NA, length(buffer.sizes)),
    avg.or.mpt = rep(NA, length(buffer.sizes)),
    avg.or.10p = rep(NA, length(buffer.sizes)),
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
    
    
    ## jackknife - leave-one-out cross-validation
    # spatial blocking by specified range with random assignment
    block <- get.jackknife(pres, back)
    
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
    testTable$or.mpt <- NA
    testTable$or.10p <- NA
    testTable$fold.auc <- NA
    
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
                                               args=c("outputformat=raw")) 
      
      test.pred <- predict(m, 
                          testSet, 
                          args=c("outputformat=cloglog")) 
      
      testTable$pred.cloglog[test.rows] <- test.pred
      
      #omission rate
      train.pred <- predict(m, 
                            trainSet, 
                            args=c("outputformat=cloglog")) 
     
      
      mpt <- min(train.pred[trainSet[,7]==1])
      p10 <- quantile(train.pred[trainSet[,7]==1], 0.1) 
      
      pred.auc <- prediction(
        c(test.pred, train.pred),
        c(1, rep(0, nrow(trainSet)))
      )
      
      testTable$auc[test.rows] <- performance(pred.auc,'auc')@y.values[[1]]
      testTable$or.mpt[test.rows] <- ifelse(test.pred < mpt, 1, 0) 
      testTable$or.10p[test.rows] <- ifelse(test.pred < p10, 1, 0) 
      
    }
    
    eval.df$avg.auc[b] <- mean(testTable$auc, na.rm = T)
    eval.df$avg.or.mpt[b] <- mean(testTable$or.mpt, na.rm = T)
    eval.df$avg.or.10p[b] <- mean(testTable$or.10p, na.rm = T)
    eval.df$npres[b] <- sum(mydata$Species==1)
    eval.df$buff.sz[b] <- buffer.sizes[b]/1000
    eval.df$species[b] <- spp.to.exerc[i]
    
    
    
  }
  eval.list[[i]] <- eval.df
  cat(paste(i, "of", length(spp.to.exerc)), fill = T)
}

eval.df.auc.1 <- eval.df
evaluation.dataset <- bind_rows(eval.list)

summ.result <- 
  evaluation.dataset %>% 
  group_by(buff.sz) %>% 
  summarise(auc.mean = mean(avg.auc, na.rm = T),
            auc.se = sqrt(sd(avg.auc, na.rm = T)),
            or.mpt.mean = mean(avg.or.mpt, na.rm = T),
            or.mpt.se = sqrt(sd(avg.or.mpt, na.rm = T)),
            or.10p.mean = mean(avg.or.10p, na.rm = T),
            or.10p.se = sqrt(sd(avg.or.10p, na.rm = T)),
            ) %>% 
  tidyr::drop_na()

(eval.res.plot<-
summ.result %>% 
  tidyr::pivot_longer(auc.mean:or.10p.se, names_to = "metric") %>% 
  filter(metric %in% c("or.mpt.mean","or.10p.mean" )) %>% 
  ggplot() +
  geom_path(aes(x = buff.sz, y = value, color = metric)) +
  geom_point(aes(x = buff.sz, y = value, color = metric)) +
  scale_color_discrete(
    name = "Threshold for omission rate", 
    labels = c("10 percentile", "Minimum presence")) +
  #geom_line(aes(x = buff.sz, y = or.mpt.mean), size = 0.2) +
  scale_x_continuous(breaks = seq(100, 1000, by = 100)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "Buffer size (Km)", y = "Omission rate") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2)
  )
)

ggsave(
  here::here("output/fig/01_buffer_size_eval_LOO.png"),
  eval.res.plot,
  height = 5,
  width = 7,
  units = "in", 
  device = "png"
)
