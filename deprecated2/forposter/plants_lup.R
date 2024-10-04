#' This script plot the symbiosis traits

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))

set.seed(1)
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))

lupulinas <- plants %>%
    filter(exp_plant == "lupulina") %>%
    filter(exp_id != "control") %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    drop_na(shoot_biomass_mg)

sam <- lupulinas %>%
    group_by(population, site_group) %>%
    count()


p <- lupulinas %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_biomass_mg, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_biomass_mg, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    geom_text(data = sam, aes(x = site_group, label = paste0("N=", n)), y = Inf, vjust = 1.2, size = 3) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    scale_y_continuous(expand = c(.1,.1)) +
    facet_wrap(.~population, scales = "free") +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(color = "grey10", fill = NA),
        strip.text = element_blank(),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot biomass (mg)")

#p <- plot_grid(p1, p2, nrow = 1, align = T, axis = "tb", scale = 0.95, rel_widths = c(1, 1))
ggsave(here::here("forposter/plants_lup.pdf"), p, width = 4, height = 4)

