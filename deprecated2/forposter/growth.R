#' This script plot the growth traits

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))

set.seed(1)
p <- gts %>%
    left_join(isolates) %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = r_30c, fill = site_group)) +
    geom_jitter(aes(x = site_group, y = r_30c, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    facet_grid(.~population, scales = "free_x") +
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
    labs(y = "growth rate at 30C (1/hr)")

#p <- plot_grid(p1, p2, nrow = 1, align = T, axis = "tb", scale = 0.95, rel_widths = c(1, 1))
ggsave(here::here("forposter/growth.pdf"), p, width = 4, height = 4)

