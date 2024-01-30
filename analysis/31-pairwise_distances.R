#' This script joins the phenotypic distance

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

dist_genetics <- read_csv(paste0(folder_data, 'temp/19-dist_genetics.csv'))
dist_traits <- read_csv(paste0(folder_data, "temp/29-dist_traits.csv"))

# Join the distances in genetics and traits
dists <- dist_genetics %>% left_join(dist_traits)
dists_long <- dists %>%
    pivot_longer(cols = c(-genome_id1, -genome_id2), names_to = "d_type", names_prefix = "d_")


write_csv(dists, paste0(folder_data, "temp/31-dists.csv"))
write_csv(dists_long, paste0(folder_data, "temp/31-dists_long.csv"))
