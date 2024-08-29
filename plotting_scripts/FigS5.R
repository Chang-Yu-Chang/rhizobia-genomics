#' This script plots the nonsymbiotic strains

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))

set.seed(1)
background_df <- tibble(exp_plant = rep(c("sativa", "lupulina"), each = 3), symbiovar = rep(c("control", "symbiotic", "nonsymbiotic"), 2))

#
plants_iso <- plants %>%
    left_join(select(iso, genome_id = genome, contig_species)) %>%
    mutate(genome_id = factor(genome_id, c(isolates$genome_id, "control"))) %>%
    mutate(symbiovar = case_when(
        contig_species %in% c("E. meliloti", "E. medicae") ~ "symbiotic",
        genome_id == "control" ~ "control",
        T ~ "nonsymbiotic"
    )) %>%
    drop_na(genome_id)

# lupulina nodule number
p1 <- plants_iso %>%
    filter(exp_nitrogen == "without nitrogen") %>%
    filter(exp_plant == "lupulina") %>%
    ggplot() +
    geom_rect(data = filter(background_df, exp_plant == "lupulina"), aes(fill = exp_plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_text(data = filter(background_df, exp_plant == "lupulina", symbiovar == "symbiotic"), aes(label = paste("M.", exp_plant)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, fontface = "italic") +
    geom_boxplot(aes(x = genome_id, y = nodule_number), outlier.size = -1, fill = NA) +
    geom_jitter(aes(x = genome_id, y = nodule_number), width = 0.05, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(.~symbiovar, scale = "free", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.text.x = element_text(angle = 15, size = 10, hjust = 0, vjust = 0 ),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "genome", y = "nodule number")


# lupulina biomass
p2 <- plants_iso %>%
    filter(exp_nitrogen == "without nitrogen") %>%
    filter(exp_plant == "lupulina") %>%
    ggplot() +
    geom_rect(data = filter(background_df, exp_plant == "lupulina"), aes(fill = exp_plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_text(data = filter(background_df, exp_plant == "lupulina", symbiovar == "symbiotic"), aes(label = paste("M.", exp_plant)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, fontface = "italic") +
    geom_boxplot(aes(x = genome_id, y = shoot_biomass_mg), outlier.size = -1, fill = NA) +
    geom_jitter(aes(x = genome_id, y = shoot_biomass_mg), width = 0.05, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(.~symbiovar, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "genome", y = "shoot biomass (mg)")

# sativa nodule number
p3 <- plants_iso %>%
    filter(exp_nitrogen == "without nitrogen") %>%
    filter(exp_plant == "sativa") %>%
    ggplot() +
    geom_rect(data = filter(background_df, exp_plant == "sativa"), aes(fill = exp_plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_text(data = filter(background_df, exp_plant == "sativa", symbiovar == "symbiotic"), aes(label = paste("M.", exp_plant)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, fontface = "italic") +
    geom_boxplot(aes(x = genome_id, y = nodule_number), outlier.size = -1, fill = NA) +
    geom_jitter(aes(x = genome_id, y = nodule_number), width = 0.05, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(.~symbiovar, scale = "free", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.text.x = element_text(angle = 15, size = 10, hjust = 0, vjust = 0 ),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "genome", y = "nodule number")

# sativa height
p4 <- plants_iso %>%
    filter(exp_nitrogen == "without nitrogen") %>%
    filter(exp_plant == "sativa") %>%
    ggplot() +
    geom_rect(data = filter(background_df, exp_plant == "sativa"), aes(fill = exp_plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_text(data = filter(background_df, exp_plant == "sativa", symbiovar == "symbiotic"), aes(label = paste("M.", exp_plant)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, fontface = "italic") +
    geom_boxplot(aes(x = genome_id, y = shoot_height), outlier.size = -1, fill = NA) +
    geom_jitter(aes(x = genome_id, y = shoot_height), width = 0.05, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(.~symbiovar, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off"
    ) +
    guides(fill = "none") +
    labs(x = "genome", y = "shoot biomass (mm)")

p <- plot_grid(p1, p2, p3, p4, ncol = 1, axis = "lr", align = "v", scale = .95, labels = c("A", "", "B", ""), rel_heights = c(1.2,1,1.2,1)) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/FigS5.png"), p, width = 8, height = 8)
