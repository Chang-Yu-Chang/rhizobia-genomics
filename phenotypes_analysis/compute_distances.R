#' This script computes the pairwise distance in traits

renv::load()
library(tidyverse)
library(janitor)
library(geosphere) # For computing the geographic distance between two locations
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv"))
gts <- read_csv(paste0(folder_phenotypes, "growth/gts.csv"))
plants_long <- read_csv(paste0(folder_phenotypes, "symbiosis/plants_long.csv"))

# 0. clean up; each row is an isolates
# 0.1. growth traits
isolates_gt <- gts %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD))

# 0.2. symbiosis traits
isolates_st <- plants_long %>%
    group_by(exp_id, trait) %>%
    summarize(value = mean(value, na.rm = T)) %>%
    pivot_wider(id_cols = exp_id, names_from = trait, values_from = value) %>%
    rename(exp_id = exp_id) %>%
    ungroup()


# 0.3 join the data
isolates_traits <- isolates %>%
    select(exp_id, genome_id, site) %>%
    left_join(select(sites, site, latitude_dec, longitude_dec)) %>%
    left_join(isolates_gt) %>%
    left_join(isolates_st)

nrow(isolates_traits) # 32 strains with growth rate data
write_csv(isolates_traits, paste0(folder_phenotypes, "isolates_traits.csv"))

# 1. Calculate pairwise euclidean distance
temp_mapping <- tibble(genome_id = isolates_traits$genome_id, item = factor(1:nrow(isolates_traits)))
compute_dist <- function (isolates_t, d_name) {
    #' This function computes the pairwise distance
    #' Each column should be one trait data
    #' Each row is one sample
    dist_t <- isolates_t %>%
        # Scale the traits
        scale() %>%
        dist(method = "euclidean", diag = T) %>%
        broom::tidy() %>%
        left_join(rename(temp_mapping, genome_id1 = genome_id, item1 = item)) %>%
        left_join(rename(temp_mapping, genome_id2 = genome_id, item2 = item)) %>%
        select(genome_id1, genome_id2, distance) %>%
        # Normalize
        mutate(distance = distance / max(distance, na.rm = T))

    colnames(dist_t)[3] <- d_name
    return(dist_t)
}

# For the 32 strains with only growth data
dist_gc <- isolates_traits %>%
    select(starts_with(c("r_", "lag_", "maxOD"))) %>%
    compute_dist("d_growth")

# For the 12 strains with both symbiosis and growth traits
dist_sym <- isolates_traits %>%
    select(shoot_biomass_mg, root_biomass_mg, nodule_count) %>%
    compute_dist("d_symbiosis") %>%
    drop_na()

# For all individual traits
list_traits <- c("shoot_biomass_mg", "root_biomass_mg", "nodule_count",
        str_subset(colnames(isolates_traits), "r_"),
        str_subset(colnames(isolates_traits), "lag_"),
        str_subset(colnames(isolates_traits), "maxOD_"))

list_dists <- rep(list(NA), length(list_traits))
for (i in 1:length(list_traits)) list_dists[[i]] <- isolates_traits[,list_traits[i]] %>% compute_dist(paste0("d_", list_traits[i]))

dist_traits <- dist_gc %>% left_join(dist_sym)
for (i in 1:length(list_traits)) dist_traits <- left_join(dist_traits, list_dists[[i]])


# 2. calculate distance by geography
dist_traits <- dist_traits %>%
    rowwise() %>%
    mutate(d_geo = distVincentySphere(
        isolates_site[match(genome_id1, isolates_site$genome_id), c("latitude_dec", "longitude_dec")],
        isolates_site[match(genome_id2, isolates_site$genome_id), c("latitude_dec", "longitude_dec")]
        ) / 1000 # convert from m to km
    ) %>%
    ungroup() %>%
    select(genome_id1, genome_id2, d_growth, d_symbiosis, d_geo, everything()) %>%
    filter(genome_id1 != genome_id2)

nrow(dist_traits) # choose(32,2) = 496
write_csv(dist_traits, paste0(folder_phenotypes, "dist_traits.csv"))
