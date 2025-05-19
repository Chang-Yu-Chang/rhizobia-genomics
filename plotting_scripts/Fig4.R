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
    mutate(contig_species = ifelse(is.na(contig_species), "control", contig_species)) %>%
    filter(population == "VA" | exp_id == "cyc") %>%
    drop_na(shoot_biomass_g) %>%
    mutate(total_biomass = shoot_biomass_g + root_biomass_g)
tb_bg <- distinct(plants, contig_species, genome_id)


# Panel A. shoot biomass----
plants_mean <- plants %>%
    group_by(contig_species, genome_id) %>%
    summarize(
        mean_biomass = mean(total_biomass), sem_biomass = sd(total_biomass) / sqrt(n()),
        mean_shoot = mean(shoot_biomass_g), sem_shoot = sd(shoot_biomass_g) / sqrt(n()),
        mean_root = mean(root_biomass_g), sem_root = sd(root_biomass_g) / sqrt(n()),
        mean_nodules = mean(nodules), sem_nodules = sd(nodules) / sqrt(n())
    )


p1 <- plants %>%
    ggplot() +
    geom_rect(data = tb_bg, aes(fill = contig_species), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .1) +
    geom_jitter(aes(x = genome_id, y = nodules), width = .1, height = 0, alpha = .3, shape = 16) +
    geom_point(data = plants_mean, aes(x = genome_id, y = mean_nodules)) +
    geom_errorbar(data = plants_mean, aes(x = genome_id, ymin = mean_nodules - qnorm(0.975) * sem_nodules, ymax = mean_nodules + qnorm(0.975) * sem_nodules), width = .1) +
    scale_fill_manual(values = species_colors) +
    facet_nested(
        contig_species+genome_id~., scales = "free_y", space = "free_y", nest_line = element_line(color = "grey90", linewidth = 2),
        strip = strip_nested(text_y = element_text(angle = 0, hjust = 0.5), background_y = element_blank()),
        switch = "y"
    ) +
    coord_flip(clip = "off") +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.clip = "off",
        strip.placement = "outside",
        panel.grid.major.x = element_line(color = "gray95"),
        panel.spacing.y = unit(0, "mm"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = unit(c(0,0,0,5), "mm")
    ) +
    guides(fill = "none") +
    labs(x = "", y = "Num. of nodules")

p2 <- plants %>%
    ggplot() +
    geom_rect(data = tb_bg, aes(fill = contig_species), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .1) +
    geom_jitter(aes(x = genome_id, y = total_biomass), width = .1, height = 0, alpha = .3, shape = 16) +
    geom_point(data = plants_mean, aes(x = genome_id, y = mean_biomass)) +
    geom_errorbar(data = plants_mean, aes(x = genome_id, ymin = mean_biomass - qnorm(0.975) * sem_biomass, ymax = mean_biomass + qnorm(0.975) * sem_biomass), width = .1) +
    scale_fill_manual(values = species_colors) +
    facet_nested(contig_species+genome_id~., scales = "free_y", space = "free_y") +
    coord_flip(clip = "off") +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray95"),
        panel.spacing.y = unit(0, "mm"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = unit(c(0,0,0,0), "mm"),
        strip.text = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "Strain", y = "Biomass (g)")

p <- plot_grid(
    p1, p2, nrow = 1,
    labels = LETTERS[1:2], rel_widths = c(2, 1),
    scale = .9,
    align = "h", axis = "tb"
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig4.png"), p, width = 6, height = 3)
