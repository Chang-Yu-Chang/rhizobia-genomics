#

library(tidyverse)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))
set.seed(42)

iso <- read_csv(paste0(folder_data, "output/iso.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv")) %>%
    left_join(select(iso, exp_id, contig_species)) %>%
    mutate(exp_id = case_when(exp_id == "control"~str_sub(exp_waterblock, 1,3), T ~ exp_id)) %>%
    filter(exp_plant == "lupulina") %>%
    mutate(contig_species = ifelse(is.na(contig_species), "control", contig_species))


#
p1 <- plants %>%
    ggplot(aes(x = exp_id, y = shoot_biomass_g)) +
    geom_boxplot(aes(fill = population), outlier.size = -1) +
    #geom_jitter(width = .1, height = 0) +
    facet_grid(~contig_species, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides() +
    labs()

p2 <- plants %>%
    #select(population, contig_species, exp_id, shoot_biomass_g) %>%
    ggplot(aes(x = exp_id, y = nodules)) +
    geom_boxplot(aes(fill = population)) +
    #geom_jitter(width = .1, height = 0) +
    facet_grid(~contig_species, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides() +
    labs()

p <- plot_grid(p1, p2, ncol = 1)
ggsave(here::here("plots/Fig4_check1.png"), p, width = 10, height = 6)


p1 <- plants %>%
    ggplot(aes(x = exp_id, y = shoot_biomass_g)) +
    geom_boxplot(aes(group = exp_waterblock, fill = population), outlier.size = -1) +
    #geom_jitter(width = .1, height = 0) +
    facet_grid(~contig_species, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides() +
    labs()

p2 <- plants %>%
    #select(population, contig_species, exp_id, shoot_biomass_g) %>%
    ggplot(aes(x = exp_id, y = nodules)) +
    geom_boxplot(aes(group = exp_waterblock, fill = population)) +
    #geom_jitter(width = .1, height = 0) +
    facet_grid(~contig_species, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides() +
    labs()

p <- plot_grid(p1, p2, ncol = 1)
ggsave(here::here("plots/Fig4_check2.png"), p, width = 10, height = 6)
