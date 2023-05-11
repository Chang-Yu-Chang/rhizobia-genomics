#' This script analyses the hand measured phenotypes: dry weight, nodule count, root weight, nodule weight

library(tidyverse)
library(broom)
library(janitor)
source(here::here("analysis/00-metadata.R"))

#
treatments <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    rename(dry_weight = `DryWeight (mg)`) %>%
    clean_names()

features <- read_csv(paste0(folder_data, "raw/rhizobia/05-root_architecture/features.csv"), show_col_types = F) %>%
    clean_names() %>%
    mutate(id = str_replace(file_name, ".png", "") %>% as.numeric())


treatments <- treatments %>% left_join(features)

#write_csv(treatments, paste0(folder_data, "temp/03-treatments.csv"))

#
treatments <- read_csv(paste0(folder_data, "temp/03-treatments.csv"), show_col_types = F)
rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), unique(treatments$rhizobia))
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")


# 1. Traits by strain ----
traits <- c("dry_weight", "nodule_number", "number_of_root_tips", "number_of_branch_points",
            "total_root_length_px", "branching_frequency_per_px", "network_area_px2",
            "average_diameter_px", "median_diameter_px", "maximum_diameter_px",
            "perimeter_px", "volume_px3", "surface_area_px2")
names(treatments)
p <- treatments %>%
    select(-contains("range_")) %>%
    pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value") %>%
    ggplot(aes(x = rhizobia, y = value, color = rhizobia_site)) +
    geom_boxplot() +
    geom_jitter(width = .1, shape = 21) +
    scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(~trait, scales = "free_y") +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(here::here("plots/03-01-trait_strain.png"), p, width = 10, height = 8)


# 2. Traits by location ----
p <- treatments %>%
    select(-contains("range_")) %>%
    pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value") %>%
    ggplot(aes(x = rhizobia_site, y = value, color = rhizobia_site)) +
    geom_boxplot() +
    geom_jitter(width = .1, shape = 21) +
    scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(~trait, scales = "free_y") +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(here::here("plots/03-02-trait_site.png"), p, width = 10, height = 8)

# 3.








