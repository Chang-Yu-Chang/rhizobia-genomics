#' This script plots the tradeoff figure

library(tidyverse)
library(janitor)
library(broom)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

iso <- read_csv(paste0(folder_data, "output/iso.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv")) %>%
    left_join(select(iso, exp_id, contig_species)) %>%
    mutate(exp_id = case_when(exp_id == "control"~str_sub(exp_waterblock, 1,3), T ~ exp_id)) %>%
    filter(exp_plant == "lupulina") %>%
    mutate(contig_species = ifelse(is.na(contig_species), "control", contig_species)) %>%
    filter(population == "VA" | exp_id == "cyc") %>%
    drop_na(shoot_biomass_g) %>%
    mutate(total_biomass = shoot_biomass_g + root_biomass_g)
plants_mean <- plants %>%
    group_by(contig_species, genome_id) %>%
    summarize(
        mean_biomass = mean(total_biomass), sem_biomass = sd(total_biomass) / sqrt(n()),
        mean_shoot = mean(shoot_biomass_g), sem_shoot = sd(shoot_biomass_g) / sqrt(n()),
        mean_root = mean(root_biomass_g), sem_root = sd(root_biomass_g) / sqrt(n()),
        mean_nodules = mean(nodules), sem_nodules = sd(nodules) / sqrt(n())
    )

gts <- read_csv(paste0(folder_phenotypes, 'growth/gts.csv')) # Growth traits per isolate
gtsl <- gts %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(temperature, exp_id, r, lag, maxOD) %>%
    mutate(loglag = -log(lag)) %>%
    pivot_longer(-c(temperature, exp_id), names_to = "trait") %>%
    left_join(select(iso, exp_id, genome_id, contig_species)) %>%
    drop_na(value) %>%
    select(-exp_id)
gtsl_mean <- gtsl %>%
    unite(trait, c(trait, temperature), sep = "_") %>%
    pivot_wider(names_from = trait)


traits_mean <- plants_mean %>%
    left_join(gtsl_mean) %>%
    filter(genome_id != "control") %>%
    select(!starts_with("sem_")) %>%
    pivot_longer(-c(contig_species, genome_id)) %>%
    group_by(name) %>%
    # Scale
    mutate(value = (value - min(value, na.rm = T))/(max(value, na.rm = T) - min(value, na.rm = T))) %>%
    mutate(trait_type = ifelse(str_detect(name, "\\dc"), "growth", "symbiosis")) %>%
    # Clean the trait name
    filter(!str_detect(name, "^log")) %>%
    mutate(trait = str_remove(name, "mean_|_\\d+c$")) %>%
    mutate(name = str_remove(name, "r_|lag_|maxOD_|mean_")) %>%
    filter(!trait %in% c("root", "shoot")) %>%
    mutate(trait = factor(case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "maxOD" ~ "yield (O.D.[600nm])",
        T ~ trait
    ), c("growth rate (1/hr)", "lag time (hr)", "yield (O.D.[600nm])", "biomass", "nodules"))) %>%
    mutate(trait = str_remove(trait, "biomass|nodules"))

# Panel A. heatmap ----
p1 <- traits_mean %>%
    ggplot() +
    geom_tile(aes(x = name, y = genome_id, fill = value), color = "black", linewidth = .5) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient2(mid = "snow", high = "maroon", low = "steelblue", midpoint = .5, breaks = c(0,.5,1)) +
    facet_nested(
        contig_species+genome_id~trait_type+trait, scales = "free", space = "free",
        switch = "y",
        strip = strip_nested(text_y = element_text(angle = 0)),
        nest_line = element_line(color = "grey80", linewidth = 1)
    ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        panel.border = element_rect(color = "black", fill = NA),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0, 5, 0, 2), "mm")
    ) +
    guides() +
    labs(x = "", y = "")



#  ----
p <- plot_grid(
    p1,
    nrow = 1, rel_widths = c(2, 1)
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig5.png"), p, width = 7, height = 4)

