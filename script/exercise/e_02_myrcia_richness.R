library(raster)
library(here)


myrcia.good <- list.files(here("output", "models", "raster_good_models"), full.names = T)
myrcia.good.LOO <- list.files(here("output", "models", "LOO", "raster_good_models"), full.names = T)

myrcia.stk <- raster::stack(c(myrcia.good, myrcia.good.LOO))

myrcia.rich <- sum(myrcia.stk, na.rm = T)

plot(myrcia.rich)

myrcia.bin <- list.files(here("output", "models", "raster_binary_prediction"), full.names = T)
myrcia.stk <- raster::stack(myrcia.bin)

myrcia.rich <- sum(myrcia.stk, na.rm = T)

plot(myrcia.rich)

myrcia.bin05 <- list.files(here("output", "models", "raster_binary_prediction", 
                                "bin_05_degree"),
                         full.names = T)

myrcia.bin05.stk <- raster::stack(myrcia.bin05)

myrcia.bin05.rich <- sum(myrcia.bin05.stk, na.rm = T)

plot(myrcia.bin05.rich)
