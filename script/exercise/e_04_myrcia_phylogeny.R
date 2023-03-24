
# load packages -----------------------------------------------------------

library(ape)
library(ade4)
library(tidyverse)
library(here)
library(taxize)

# load data ---------------------------------------------------------------

# species names

names_distribution <- 
list.files(here("output", "models", "raster_binary_prediction", "bin_05_degree")) %>% 
  word(1, 2, "_")

# tree

myrcia_tree <- read.tree(here("data/extra/myrcia_group/myrcia_clean.txt"))

phy_names <- myrcia_tree$tip.label %>% 
  tolower()


# find species macth ------------------------------------------------------


pow.review <- data.frame(
  name_searched = character(),
  taxonomic_status = character(),
  nomenclatural_status = character(),
  name_accepted = character()
)

phy_names_2 <- str_replace_all(phy_names,"_aff_", " ") %>% 
  str_replace_all( "_", " ") %>% 
  str_to_sentence() 

for(i in 1:length(phy_names_2)){
  x <- phy_names_2[i]
  
  g <- get_pow(x, ask = T, messages = F)
  pow_id <- as.character(g)
  if(is.na(g)){
    found <- attr(g,"match")
    multi <- ifelse(attr(g,"multiple_matches"), "multiple_matches", NA)
    
    add_pow <- data.frame(
      name_searched = x,
      taxonomic_status = found,
      nomenclatural_status = multi,
      name_accepted = NA
    ) 
    
    pow.review <- add_row(
      pow.review, add_pow
    )
    
    next()
  }
  
  lk <- pow_lookup(pow_id)
  
  tax_status <- lk$meta$taxonomicStatus
  nom_status <- lk$meta$nomenclaturalStatus
  
  if(tax_status != "Accepted"){
    if(tax_status == "Unplaced"){
      accepted.name <- NA
    }else{
      accepted.name <- lk$meta$accepted$name
    }
    
  } else {
    accepted.name <- x
  }
  
  add_pow <- data.frame(
    name_searched = x,
    taxonomic_status = tax_status,
    nomenclatural_status = nom_status,
    name_accepted = accepted.name
  ) 
  
  pow.review <- add_row(
    pow.review, add_pow
  )
  
}

sum(
pow.review$name_accepted %>% 
  tolower() %>% 
  str_replace_all(" ", "_") %in% names_distribution
)

pow.review <- 
  pow.review %>% 
  mutate(
    phy_names = phy_names,
    str_order = 1:nrow(pow.review))

pow.resolved <- 
  pow.review %>% 
  filter(nomenclatural_status != "multiple_matches",
         taxonomic_status != "not found")



pow.not.found <- 
  pow.review %>% 
  filter(taxonomic_status == "not found") 

pow.not.found$name_accepted <- c(
  "Myrcia anacardiifolia", #mispelled
  NA, #dubious
  "Myrcia pendula", #mispelled
  NA, 
  "Myrcia subalpestris", #synonym
  NA, 
  NA, 
  NA,
  "Myrcia clarendonensis", #synonym
  "Myrcia plusiantha", #synonym
  NA, 
  "Myrcia megaphylla", #flora do brasil - aceito
  NA,
  "Myrcia spathulifolia", #synonym,
  rep(NA, 8)
)


pow.review.2 <- 
bind_rows(pow.resolved, pow.not.found) %>% 
  arrange(str_order) %>% 
  mutate(
    phy_new_names = ifelse(
      is.na(name_accepted), name_searched, name_accepted
      ), 
    phy_new_names = str_replace_all(phy_new_names, " ", "_") %>% 
      tolower() 
    )

write.csv(pow.review.2, here("data", "extra", "phy_names_reviewed_pow_df.csv"))
  
### clean phylogeny

## remove all tips waith taxa not identified at species level
## remove tips with duplicated revised accepted names
## - keep tip with the same name as the accepetd vesion
## - keep tip with aff to the accepted version
## - choose one tip randomly

pow.review.2$phy_new_names[pow.review.2$phy_new_names %in% names_distribution]

pow.review.2 %>% 
  filter(
    duplicated(phy_new_names) | duplicated(phy_new_names, fromLast = TRUE)
  ) %>% 
  arrange(phy_new_names)
