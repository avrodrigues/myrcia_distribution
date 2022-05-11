
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

myrcia.bin05 <- list.files(here("output", "models", "raster_binary_prediction", 
                                "bin_05_degree"),
                           full.names = T)

myrcia.stk <- stack(myrcia.bin05)

myrcia.df <- as.data.frame(myrcia.stk, xy = T)
names(myrcia.df)


myrcia.df.long <- 
myrcia.df %>% 
  as_tibble() %>% 
  mutate(site = 1:nrow(myrcia.df)) %>% 
  pivot_longer(cols = myrcia_aethusa_fc_LQP_rm_1:myrcia_zuzygium_fc_LQP_rm_2, 
               values_to = "presence", 
               names_to = "species") %>% 
  mutate(
    species = word(species, 1, 2, sep = "_")
  )

myrcia.comp.mtx <- 
  myrcia.df.long %>% 
    pivot_wider(names_from = species, values_from = presence) %>% 
    dplyr::select(!1:3) %>% 
    as.matrix()

comp3sp <- which(rowSums(myrcia.comp.mtx) >=3)

xy_comp3sp <- myrcia.df[comp3sp, 1:2]

# compute richness --------------------------------------------------------

myrcia.rich <-
myrcia.df.long %>% 
  group_by(site) %>% 
  summarise(rich = sum(presence==1))

myrcia.rich <- as.data.frame(sum(myrcia.stk), xy = T)
myrcia.rich <- myrcia.rich %>% 
  mutate(
    site = 1:nrow(myrcia.rich), 
    layer = ifelse(layer == 0, NA, layer))

# plot draft
(rich.map <- 
ggplot(myrcia.rich) +
  geom_raster(aes(x = x, y = y, fill = layer )) +
  scale_fill_gradientn(
    name = "Species Richness",
    colors=met.brewer("VanGogh3"),
    na.value = "#FAF8F4")+
  theme_bw() +
  coord_fixed()+
  theme(legend.position = "bottom")
)

ggsave(
  here("output", "fig", "richness_draft.png"),
  rich.map,
  width = 6,
  height = 6.5
)


# beta diversity turnover -------------------------------------------------


myrcia.comp.hel <- decostand(myrcia.comp.mtx[comp3sp,], "hellinger")
myrcia.dist <- dist(myrcia.comp.hel)

myrcia.pcoa <- pcoa(myrcia.dist)

axis.sel <- which(myrcia.pcoa$values$Relative_eig > 0.05)

pcoa.axes.sel <- as.data.frame(myrcia.pcoa$vectors[,axis.sel])

### regionalization (betadiversity turnover)
dist.pcoa <- dist(pcoa.axes.sel)
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

gr <- cutree(pcoa.ward, k=k.best)

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

## plot draft
beta.df <- 
bind_cols(xy_comp3sp, pcoa.axes.sel) %>% 
  mutate(k_group = gr)

my_theme <- list(
  theme(
    panel.background = element_rect(fill = "#FAF8F4"), 
    panel.grid = element_blank(), 
    text = element_text(color = greys[1]), 
    axis.text = element_text(color = greys[2]), 
    axis.ticks = element_line(color = greys[3]), 
    panel.border = element_rect(color = greys[4], fill = NA), 
    axis.title = element_blank()
  )
)


coast <- rnaturalearth::ne_coastline(returnclass = "sf")

(beta.disc <-
ggplot() +
  geom_raster(data = beta.df, aes(x = x, y = y, fill = as.factor(k_group))) +
  scale_fill_manual(values = (blue_gold_red[c(1,3,4,2)]))+
  geom_sf(data = coast, color = greys[2]) +
  coord_sf(xlim = c(-110, -35), ylim = c(-40, 25)) +
  theme(legend.position = 'none') +
  my_theme
)

ggsave(
  here("output", "fig", "beta_discrete_draft.png"),
  beta.disc,
  width = 6,
  height = 6.5
)




