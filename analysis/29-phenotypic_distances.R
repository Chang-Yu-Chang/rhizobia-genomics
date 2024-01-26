#' This script joins the phenotypic distance

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
source(here::here("analysis/00-metadata.R"))

gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'), show_col_types = F)
treatments_long <- read_csv(paste0(folder_data, "temp/23-treatments_long.csv"), show_col_types = F)

# 0. clean up; basically one sample per row ----
# 0.1. growth traits ----
isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD))

# 0.2. symbiosis traits ----
isolates_sym <- treatments_long %>%
    group_by(rhizobia, trait) %>%
    summarize(value = mean(value, na.rm = T)) %>%
    pivot_wider(id_cols = rhizobia, names_from = trait, values_from = value) %>%
    rename(exp_id = rhizobia) %>%
    ungroup()

# Clean up the rhizobia name
isolates_sym <- isolates_sym %>%
    mutate(exp_id = str_replace(exp_id, "b_", "_")) %>%
    mutate(exp_id = str_remove(exp_id, "_c\\d")) %>%
    mutate(exp_id = str_replace(exp_id, "p_", "p")) %>%
    mutate(exp_id = str_replace(exp_id, "_", "-"))

# 0.3 join the data
isolates_traits <- isolates_mapping %>%
    select(exp_id, genome_id) %>%
    left_join(isolates_gc) %>%
    left_join(isolates_sym)

# 0.4 calculate pairwise euclidean distance
# 32 strains
temp_mapping <- tibble(genome_id = isolates_traits$genome_id, item = factor(1:nrow(isolates_traits)))
dist_gc <- isolates_traits %>%
    select(starts_with(c("r_", "lag_", "maxOD"))) %>%
    # Scale the traits
    scale() %>%
    dist(method = "euclidean") %>%
    broom::tidy() %>%
    left_join(rename(temp_mapping, genome_id1 = genome_id, item1 = item)) %>%
    left_join(rename(temp_mapping, genome_id2 = genome_id, item2 = item)) %>%
    rename(distance_growth = distance) %>%
    select(genome_id1, genome_id2, distance_growth) %>%
    # Normalize
    mutate(distance_growth = distance_growth / max(distance_growth, na.rm = T))

# 12 strains
dist_sym <- isolates_traits %>%
    select(dry_weight_mg, nodule_number, root_weight_mg) %>%
    # Scale the traits
    scale() %>%
    dist(method = "euclidean") %>%
    broom::tidy() %>%
    left_join(rename(temp_mapping, genome_id1 = genome_id, item1 = item)) %>%
    left_join(rename(temp_mapping, genome_id2 = genome_id, item2 = item)) %>%
    rename(distance_symbiosis = distance) %>%
    select(genome_id1, genome_id2, distance_symbiosis) %>%
    # Normalize
    mutate(distance_symbiosis = distance_symbiosis / max(distance_symbiosis, na.rm = T)) %>%
    drop_na


# 0.5 join the distances----
dists_trait <- dist_gc %>% left_join(dist_sym)

write_csv(dists_trait, paste0(folder_data, "temp/29-dists_trait.csv"))


# 1. phenotpyic distances ----
p <- dists_trait %>%
    drop_na() %>%
    ggplot() +
    geom_point(aes(x = distance_growth, y = distance_symbiosis), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/29-01-trait_dist.png"), p, width = 4, height = 4)


#
cor.test(dists_trait$distance_growth, dists_trait$distance_symbiosis)














