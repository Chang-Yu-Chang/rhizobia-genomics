#' This script cleans the tables of phenotypic distance and join them

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'), show_col_types = F)
plants_long <- read_csv(paste0(folder_data, "temp/23-plants_long.csv"), show_col_types = F)

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

# 0.3 join the data
isolates_traits <- isolates_mapping %>%
    select(exp_id, genome_id) %>%
    left_join(isolates_gc) %>%
    left_join(isolates_sym)

nrow(isolates_traits) # 32 strains with growth rate data 

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

write_csv(dist_traits, paste0(folder_data, "temp/24-dist_traits.csv"))
