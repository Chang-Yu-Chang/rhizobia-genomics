#' This script cleans the tables of phenotypic distance and join them

renv::load()
library(tidyverse)
library(janitor)
library(geosphere) # For computing the geographic distance between two locations
source(here::here("analysis/00-metadata.R"))

isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"))
sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"))
gc_prm_summs <- read_csv(paste0(folder_data, "temp/21-gc_prm_summs.csv"))
plants_long <- read_csv(paste0(folder_data, "temp/23-plants_long.csv"))

# 0. clean up; basically one sample per row
# 0.1. growth traits
isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD))

# 0.2. symbiosis traits
isolates_sym <- plants_long %>%
    group_by(exp_id, trait) %>%
    summarize(value = mean(value, na.rm = T)) %>%
    pivot_wider(id_cols = exp_id, names_from = trait, values_from = value) %>%
    rename(exp_id = exp_id) %>%
    ungroup()

# 0.3 coordinate of sites
isolates_site <- isolates_mapping %>% left_join(sites)

# 0.3 join the data
isolates_traits <- isolates_mapping %>%
    select(exp_id, genome_id) %>%
    left_join(isolates_gc) %>%
    left_join(isolates_sym) %>%
    left_join(isolates_site)

nrow(isolates_traits) # 32 strains with growth rate data 

write_csv(isolates_traits, paste0(folder_data, "temp/24-isolates_traits.csv"))

# 1. Calculate pairwise euclidean distance
# For the 32 strains with only growth data
temp_mapping <- tibble(genome_id = isolates_traits$genome_id, item = factor(1:nrow(isolates_traits)))
dist_gc <- isolates_traits %>%
    select(starts_with(c("r_", "lag_", "maxOD"))) %>%
    # Scale the traits
    scale() %>%
    dist(method = "euclidean", diag = T) %>%
    broom::tidy() %>%
    left_join(rename(temp_mapping, genome_id1 = genome_id, item1 = item)) %>%
    left_join(rename(temp_mapping, genome_id2 = genome_id, item2 = item)) %>%
    rename(d_growth = distance) %>%
    select(genome_id1, genome_id2, d_growth) %>%
    # Normalize
    mutate(d_growth = d_growth / max(d_growth, na.rm = T))

# For the 12 strains with both symbiosis and growth traits
dist_sym <- isolates_traits %>%
    select(dry_weight_mg, nodule_number, root_weight_mg) %>%
    # Scale the traits
    scale() %>%
    dist(method = "euclidean", diag = T) %>%
    broom::tidy() %>%
    left_join(rename(temp_mapping, genome_id1 = genome_id, item1 = item)) %>%
    left_join(rename(temp_mapping, genome_id2 = genome_id, item2 = item)) %>%
    rename(d_symbiosis = distance) %>%
    select(genome_id1, genome_id2, d_symbiosis) %>%
    # Normalize
    mutate(d_symbiosis = d_symbiosis / max(d_symbiosis, na.rm = T)) %>%
    drop_na

# 
dist_traits <- dist_gc %>% left_join(dist_sym)

# 2. calculate distance by geography
dist_traits <- dist_traits %>%
    rowwise() %>%
    mutate(d_geo = distVincentySphere(
        isolates_site[match(genome_id1, isolates_site$genome_id), c("latitude_dec", "longitude_dec")],
        isolates_site[match(genome_id2, isolates_site$genome_id), c("latitude_dec", "longitude_dec")]
        ) / 1000 # convert from m to km
    )

write_csv(dist_traits, paste0(folder_data, "temp/24-dist_traits.csv"))


