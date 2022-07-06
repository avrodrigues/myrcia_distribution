
# load packages -----------------------------------------------------------

library(raster)
library(tidyverse)
library(ape)
library(vegan)
library(here)
library(MetBrewer)
library(cluster)
library(rnaturalearth)


# load data ---------------------------------------------------------------

l.myrcia.df <- 
  readRDS(here("output", "models", "list_myrcia_binary_df_05_degree.rds"))

names(l.myrcia.df)

l.myrcia.wide <- 
  map(
    l.myrcia.df, 
    function(myrcia.df){
      
      myrcia.df %>% 
        pivot_wider(
          names_from = species, 
          values_from = presence, 
          values_fill = 0) 
      
    }
  )

l.myrcia.comp.mtx <- 
  map(
    l.myrcia.wide, 
    function(myrcia.wide){
      myrcia.wide %>% 
        dplyr::select(!1:2) %>% 
        as.matrix()
    }
  )
  

l_comp3sp <- map(l.myrcia.comp.mtx, function(x){
  which(rowSums(x) >=3)
})
  
  

l_xy_comp3sp <- map(1:3, function(i){
  
  xy <- 
    l.myrcia.wide[[i]] %>% 
      dplyr::select(1:2) %>% 
      as.matrix()

  comp3sp <- l_comp3sp[[i]]
    
  xy[comp3sp, ]
} )

# compute richness --------------------------------------------------------

l.myrcia.rich <- map(
  l.myrcia.df, 
  function(myrcia.df) {
    myrcia.df %>% 
      group_by(x, y) %>% 
      summarise(rich = sum(presence))
  }
)


# |- map theme ------------------------------------------------------------

five_hues <- c(
  "#3d291a",
  "#a9344f",
  "#578a5b",
  "#83a6c4",
  "#fcc573"
)

blue_gold_red <- c(
  "#303260",
  "#88beca",
  "#e4b434",
  "#7b5a28",
  "#e19296",
  "#c62f22"
  
)

greys <- c(
  "#040400",
  "#1F1F1B",
  "#3B3B37",
  "#575753",
  "#73736F",
  "#8E8E8A",
  "#AAAAA6",
  "#C6C6C2",
  "#E2E2DE",
  "#FEFEFA"
)

my_theme <- list(
  theme(
    panel.background = element_rect(fill = "#FAF8F4"), 
    panel.grid = element_blank(), 
    text = element_text(color = greys[1]), 
    title = element_text(color = greys[2]),
    axis.text = element_text(color = greys[2]), 
    axis.ticks = element_line(color = greys[3]), 
    panel.border = element_rect(color = greys[4], fill = NA), 
    axis.title = element_blank()
  )
)
coast <- rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)



# |- map richness --------------------------------------------------------

max_rich <- map_dbl(l.myrcia.rich, function(x){
  max(x$rich)
}) %>% max()

l.rich.map <- map(l.myrcia.rich, function(myrcia.rich){
  ggplot(myrcia.rich) +
    geom_raster(aes(x = x, y = y, fill = rich )) +
    scale_fill_gradientn(
      name = "Species Richness",
      colors=met.brewer("VanGogh3"),
      na.value = "#FAF8F4", 
      limits = c(0, max_rich)
    ) +
    geom_sf(data = coast, color = greys[6]) +
    coord_sf(xlim = c(-110, -35), ylim = c(-40, 25)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    my_theme
})

l.rich.map$few_occ_thr_min +
  ggtitle("Threshold minimum presence")
l.rich.map$few_occ_thr_10p +
  ggtitle("Threshold 10 percentile presence")
l.rich.map$few_occ_thr_site +
  ggtitle("Threshold dependent on site richness")

ggsave(
  here("output", "fig", "richness_draft.png"),
  rich.map,
  width = 6,
  height = 6.5
)


# beta diversity turnover -------------------------------------------------

pcoa.axes.sel <- map(
  1:3, 
  function(i){
    comp <- l.myrcia.comp.mtx[[i]]
    keep.rows <- l_comp3sp[[i]]
    
    myrcia.dist <- vegan::vegdist(comp[keep.rows,])
    myrcia.pcoa <- pcoa(myrcia.dist)
    axis.sel <- which(myrcia.pcoa$values$Relative_eig > 0.05)
    
    as.data.frame(myrcia.pcoa$vectors[,axis.sel])
  })


### regionalization (betadiversity turnover)

l.gr <- map(pcoa.axes.sel, function(x){
  dist.pcoa <- dist(x)
  pcoa.ward <- hclust(dist.pcoa, method="average")
  
  
  pcoa.ward.coph <- cophenetic(pcoa.ward)
  cor(dist.pcoa, pcoa.ward.coph)
  
  ## Escolha do número de grupos pelo metodo de Largura Média da Silhueta
  ## 
  
  asw <- numeric(9)
  for (k in 2:10) {
    sil <- silhouette(cutree(pcoa.ward, k=k), dist.pcoa)
    asw[k-1] <- summary(sil)$avg.width
  }
  k.best <- which.max(asw)
  
 cutree(pcoa.ward, k=k.best)
})




## plot draft
l.beta.df <- map(1:3, function(i){
  bind_cols(l_xy_comp3sp[[i]], pcoa.axes.sel[[i]]) %>% 
    mutate(k_group = l.gr[[i]])
})



l.beta.map <- map(l.beta.df, function(beta.df){
  ggplot() +
  geom_raster(data = beta.df, aes(x = x, y = y, fill = as.factor(k_group))) +
  scale_fill_manual(values = (blue_gold_red[c(1,3,4,2,5)]))+
  geom_sf(data = coast, color = greys[2]) +
  coord_sf(xlim = c(-110, -35), ylim = c(-40, 25)) +
  theme(legend.position = 'none') +
  my_theme

})

l.beta.map[[3]]

ggsave(
  here("output", "fig", "beta_discrete_draft.png"),
  beta.disc,
  width = 6,
  height = 6.5
)




