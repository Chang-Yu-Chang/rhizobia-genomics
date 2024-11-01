#' This script plots the individaul plant data of plant experiments

library(tidyverse)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

# Summary
nrow(plants) # 822 plants
table(plants$exp_plant) # 251 lupulina, 571 sativa
length(unique(plants$exp_id)) # 26 rhizobia strains used in the plant experiments


# 1. All plant traits ----
set.seed(1)
p <- plants %>%
    filter(population != "control", exp_plant == "sativa", exp_nitrogen == "without nitrogen") %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    ggplot(aes(x = population, y = value, fill = population)) +
    geom_boxplot(outlier.shape = -1, alpha = 0.5) +
    geom_jitter(size = .5, shape = 21, width = .2) +
    scale_fill_manual(values = population_colors) +
    coord_cartesian(clip = "off") +
    facet_nested(
        gradient~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none",
        strip = strip_nested(bleed=T, clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("grey90", 4), rep("white", 10))))) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = "grey90"),
    ) +
    guides(fill = "none") +
    labs(y = "")

ggsave(paste0(folder_phenotypes, "plants/01-sativa_traits.png"), p, width = 16, height = 8)

# 2. Strain level data----
set.seed(1)
idf <- arrange(plants, population) %>% distinct(exp_id) %>% pull(exp_id)

p <- plants %>%
    filter(population != "control", exp_plant == "sativa", exp_nitrogen == "without nitrogen") %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    mutate(exp_id = factor(exp_id, idf)) %>%
    drop_na(exp_id) %>%
    ggplot(aes(x = exp_id, y = value, fill = population)) +
    geom_boxplot(outlier.shape = -1, alpha = 0.5) +
    geom_jitter(size = .5, shape = 21, width = .2) +
    scale_fill_manual(values = population_colors) +
    coord_cartesian(clip = "off") +
    facet_nested(
        gradient~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none",
        strip = strip_nested(bleed=T, clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("grey90", 4), rep("white", 10))))) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = "grey90"),
    ) +
    guides(fill = "none") +
    labs(y = "")

ggsave(paste0(folder_phenotypes, "plants/02-sativa_traits_strain.png"), p, width = 16, height = 8)


