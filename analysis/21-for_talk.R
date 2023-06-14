#' This script make figures for talks

library(tidyverse)
library(cowplot)
library(broom)
library(janitor)
library(ggsci)
library(waffle) #remotes::install_github("hrbrmstr/waffle")
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)

# 1. map for sampling site ----
world_coordinates <- map_data("world")

# create world map using ggplot() function
ggplot() +

    # geom_map() function takes world coordinates
    # as input to plot world map
    geom_map(
        data = world_coordinates, map = world_coordinates,
        aes(long, lat, map_id = region)
    )

# 2. factorial design ----
treatments %>% tabyl(rhizobia_site, plant_site, show_missing_levels = T)
treatments %>% tabyl(rhizobia, plant_site, show_missing_levels = T)

p <- treatments %>%
    mutate(rhizobia_site = ifelse(is.na(rhizobia_site), "control", rhizobia_site),
           rhizobia = ifelse(is.na(rhizobia), "control", rhizobia)) %>%
    #filter(rhizobia != "control") %>%
    mutate(rhizobia_site = factor(rhizobia_site)) %>%
    group_by(rhizobia_site, rhizobia, plant_site) %>%
    count(.drop = F) %>%
    arrange(rhizobia_site, plant_site) %>%
    ggplot() +
    geom_waffle(aes(fill = plant_site, values = n),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = TRUE) +
    scale_fill_manual(values = rhizobia_site_colors) +
    #scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Set1")) +
    #scale_y_continuous(labels = function(x) x * 10, expand = c(0,0)) +
    scale_alpha_manual(values = rhizobia_alphas) +
    facet_grid(~rhizobia_site) +
    coord_equal() +
    theme_void() +
    theme(
        legend.position = "bottom",
        strip.text.x = element_text(hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/21-02-factorial_design.png"), p, width = 5, height = 2)


# 2a. factorial design, by plant and rhizobia ----
cc <- c(H = "#0C6291", L = "#BF4342", control = "#CBD4C2")
p <- treatments %>%
    mutate(rhizobia_site = ifelse(is.na(rhizobia_site), "control", rhizobia_site),
           rhizobia = ifelse(is.na(rhizobia), "control", rhizobia)) %>%
    mutate(rhizobia_site = factor(rhizobia_site, c("H", "L", "control"))) %>%
    group_by(rhizobia_site, rhizobia, plant_site) %>%
    count(.drop = F) %>%
    arrange(rhizobia_site, plant_site) %>%
    ggplot() +
    geom_waffle(aes(values = n, fill = rhizobia_site),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = F) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = cc, labels = c("high", "low", "control"), breaks = c("H", "L", "control")) +
    coord_equal() +
    theme_void() +
    theme(
        legend.position = "bottom",
        strip.text.x = element_text(hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = guide_legend(title = "rhizobia site")) +
    labs()
ggsave(paste0(folder_data, "temp/21-02a-factorial_design_plain.png"), p, width = 3.5, height = 2)

# 2b. ----
cc <- c(H = "#0C6291", L = "#BF4342", control = "#CBD4C2")
p <- treatments %>%
    mutate(rhizobia_site = ifelse(is.na(rhizobia_site), "control", rhizobia_site),
           rhizobia = ifelse(is.na(rhizobia), "control", rhizobia)) %>%
    mutate(rhizobia_site = factor(rhizobia_site, c("H", "L", "control"))) %>%
    group_by(rhizobia_site, rhizobia, plant_site) %>%
    count(.drop = F) %>%
    arrange(rhizobia_site, rhizobia, plant_site) %>%
    ggplot() +
    geom_waffle(aes(values = n, fill = rhizobia_site, alpha = rhizobia),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = F) +
    #scale_alpha_manual(values = rhizobia_alphas[1:6]) +
    scale_fill_manual(values = cc, labels = c("high", "low", "control"), breaks = c("H", "L", "control")) +
    facet_wrap(~rhizobia_site) +
    coord_equal() +
    theme_void() +
    theme(
        legend.position = "bottom",
        strip.text.x = element_text(hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = guide_legend(title = "rhizobia site")) +
    labs()
ggsave(paste0(folder_data, "temp/21-02b-factorial_design_strains.png"), p, width = 3.5, height = 2)

# 3. conceptual plot ----

set.seed(1)
tb1 <- tibble(site = rep(c("site 1","site 2"), each = 20),
       value = rnorm(40))
tb2 <- tibble(site = rep(c("site 1","site 2"), each = 20),
              value = c(rnorm(20, 0), rnorm(20, 3)))


p1 <- tb1 %>%
    ggplot(aes(x = site, y = value, group = site)) +
    geom_jitter(width = 0.1, shape = 21, size = 2, stroke = 1) +
    theme_classic() +
    theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    ) +
    guides() +
    labs(x = "", y = "trait")

p2 <- tb2 %>%
    ggplot(aes(x = site, y = value, group = site)) +
    geom_jitter(width = 0.1, shape = 21, size = 2, stroke = 1) +
    theme_classic() +
    theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    ) +
    guides() +
    labs(x = "", y = "trait")

p <- plot_grid(p1, NULL, p2, nrow = 1, rel_widths = c(1,.4,1), axis = "tb", align = "h") + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/21-03-conceptual.png"), p, width = 5, height = 2.5)
