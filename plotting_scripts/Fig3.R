#' This script plots the growth traits

library(tidyverse)
library(cowplot)
library(ggh4x) # for nested facets
source(here::here("metadata.R"))

# Prepare the data ----
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    mutate(contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")))
gc_summs <- read_csv(paste0(folder_data, "phenotypes/growth/gc_summs.csv"))
gts <- read_csv(paste0(folder_phenotypes, 'growth/gts.csv')) # Growth traits per isolate
gc_summs <- left_join(gc_summs, select(iso, exp_id, genome_id, contig_species))

gtsl <- gts %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(temperature, exp_id, r, lag, maxOD) %>%
    pivot_longer(-c(temperature, exp_id), names_to = "trait") %>%
    left_join( select(iso, exp_id, genome_id, contig_species)) %>%
    drop_na(value)

# Panel A
p1 <- ggdraw()

# Panel growth curve
p2 <- gc_summs %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs, group = exp_id, color = contig_species), linewidth = .3) +
    geom_text(aes(label = temperature), x = 5, y = .45) +
    scale_color_manual(values = species_colors) +
    scale_x_continuous(breaks = seq(0, 48, 12)) +
    coord_cartesian(clip = "off") +
    facet_wrap2(~temperature, nrow = 1) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey95", linewidth = .5),
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 10),
        strip.placement = "outside",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(5, "mm"),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.background = element_blank()
    ) +
    guides() +
    labs(x = "Time (hour)", y = "O.D.[600nm]")

# Panel C
gtwlm <- gtsl %>%
    group_by(temperature, contig_species, trait) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value))) %>%
    group_by(temperature, trait) %>%
    mutate(max_mean_value = max(mean_value, na.rm = T)) %>%
    replace_na(list(ci_value = 0)) %>%
    mutate(trait = case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "maxOD" ~ "yield [OD]"
    ))

p3 <- gtsl %>%
    mutate(trait = case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "maxOD" ~ "yield [OD]"
    )) %>%
    ggplot() +
    # Each strain
    geom_line(aes(x = temperature, y = value, group = exp_id, color = contig_species), alpha = .1) +
    # Mean value
    geom_ribbon(data = gtwlm, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = contig_species, group = contig_species), inherit.aes = FALSE, alpha = 0.2) +
    geom_line(data = gtwlm, aes(x = temperature, y = mean_value, color = contig_species, group = contig_species)) +
    geom_point(data = gtwlm, aes(x = temperature, y = mean_value, color = contig_species, group = contig_species)) +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    facet_wrap2(~trait, scales = "free_y", nrow = 1, strip.position = "left") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey95", linewidth = .5),
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10),
        strip.placement = "outside",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(5, "mm"),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.background = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(color = NA), nrow = 2, direction = "vertical")) +
    labs(x = expression(paste("Temperature (", degree, "C)")))

#p_left <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,3))
p <- plot_grid(
    p2, p3,
    ncol = 1, align = "v", axis = "lr",
    labels = c("A", "B"), scale = .95
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig3.png"), p, width = 8, height = 6)
