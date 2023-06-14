#' This script plots the root traits

library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

# growth curves
gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L"))

# common garden
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
treatments_long <- read_csv(paste0(folder_data, "temp/11-treatments_long.csv"), show_col_types = F)
treatments_scaled <- read_csv(paste0(folder_data, "temp/11-treatments_scaled.csv"), show_col_types = F)
treatments_scaled_long <- read_csv(paste0(folder_data, "temp/11-treatments_scaled_long.csv"), show_col_types = F)

# 0. combine the data from treatments and gc ----
# Subset only Ensifer -----
subset_ensifer <- function(tb) {
    tb %>%
        left_join(select(isolates_RDP, strain, Genus)) %>%
        drop_na()
}

gc <- gc %>% subset_ensifer()
gc_summ <- gc_summ %>% subset_ensifer()
gc.prm <- gc.prm %>% subset_ensifer()
gc.prm.stat <- gc.prm.stat %>% subset_ensifer()

# Combine the data from rhizobia and ----
treatments_long_stat <- treatments_long %>%
    filter(trait %in% c("dry_weight_mg", "nodule_number", "root_weight_mg", "total_root_length_px", "branching_frequency_per_px", "network_area_px2", "average_diameter_px")) %>%
    group_by(rhizobia_site, rhizobia, trait) %>%
    drop_na(rhizobia, value) %>%
    summarize(trait_mean = mean(value, na.rm = T), trait_sd = sd(value, na.rm = T), n = n())

t1 <- gc.prm.stat %>%
    rename(rhizobia = strain) %>%
    filter(rhizobia %in% rhizobia_strains) %>%
    select(rhizobia, r, lag, maxOD) %>%
    pivot_longer(cols = -rhizobia, names_to = "trait", values_to = "trait_mean")
t2 <- gc.prm.stat %>%
    rename(rhizobia = strain) %>%
    filter(rhizobia %in% rhizobia_strains) %>%
    select(rhizobia, r=r.sem, lag=lag.sem, maxOD=maxOD.sem) %>%
    pivot_longer(cols = -rhizobia, names_to = "trait", values_to = "trait_sd")

gc_long_stat <- left_join(t1, t2) %>% mutate(n = 4) %>% mutate(rhizobia_site = str_sub(rhizobia, 1, 1))

trait_long_stat <- bind_rows(mutate(gc_long_stat, trait_type = "growth"), mutate(treatments_long_stat, trait_type = "mutualism")) %>%
    select(rhizobia_site, rhizobia, everything()) %>%
    arrange(rhizobia_site, rhizobia, trait)

# 1. r vs. biomass ----
p <- trait_long_stat %>%
    group_by(rhizobia_site, rhizobia) %>%
    select(rhizobia_site, rhizobia, trait, trait_mean, trait_sd) %>%
    filter(trait %in% c("dry_weight_mg", "r")) %>%
    pivot_wider(names_from = trait, values_from = c(trait_mean, trait_sd)) %>%
    ggplot() +
    geom_point(aes(x = trait_mean_r, y = trait_mean_dry_weight_mg, color = rhizobia_site, group = rhizobia), shape = 21, size = 3, stroke = 1) +
    geom_segment(aes(x = trait_mean_r, xend = trait_mean_r, y = trait_mean_dry_weight_mg + trait_sd_dry_weight_mg, yend = trait_mean_dry_weight_mg - trait_sd_dry_weight_mg, color = rhizobia_site)) +
    geom_segment(aes(x = trait_mean_r + trait_sd_r, xend = trait_mean_r - trait_sd_r, y = trait_mean_dry_weight_mg, yend = trait_mean_dry_weight_mg, color = rhizobia_site)) +
    scale_color_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/11c-01-r_vs_biomass.png"), p, width = 5, height = 3)


# 2. r vs. nodule number ----
p <- trait_long_stat %>%
    group_by(rhizobia_site, rhizobia) %>%
    select(rhizobia_site, rhizobia, trait, trait_mean, trait_sd) %>%
    filter(trait %in% c("nodule_number", "r")) %>%
    pivot_wider(names_from = trait, values_from = c(trait_mean, trait_sd)) %>%
    ggplot() +
    geom_point(aes(x = trait_mean_r, y = trait_mean_nodule_number, color = rhizobia_site, group = rhizobia), shape = 21, size = 3, stroke = 1) +
    geom_segment(aes(x = trait_mean_r, xend = trait_mean_r, y = trait_mean_nodule_number + trait_sd_nodule_number, yend = trait_mean_nodule_number - trait_sd_nodule_number, color = rhizobia_site)) +
    geom_segment(aes(x = trait_mean_r + trait_sd_r, xend = trait_mean_r - trait_sd_r, y = trait_mean_nodule_number, yend = trait_mean_nodule_number, color = rhizobia_site)) +
    scale_color_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/11c-02-r_vs_nodule.png"), p, width = 5, height = 3)













