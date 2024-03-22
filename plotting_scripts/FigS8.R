#' This script plots the reaction norm

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"))
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'))
isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_longer(cols = -c(exp_id, temperature), names_to = "trait") %>%
    left_join(isolates)

p <- isolates_gc %>%
    filter(temperature != "40c") %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    ggplot() +
    geom_point(aes(x = temperature, y = value, color = site_group, group = exp_id), shape = 21, size = 2, stroke = 1) +
    geom_line(aes(x = temperature, y = value, color = site_group, group = exp_id), alpha = 0.5) +
    scale_color_manual(values = site_group_colors) +
    facet_grid(trait ~ population, scales = "free_y") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
    ) +
    guides() +
    labs()

ggsave(here::here("plots/FigS8.png"), p, width = 8, height = 6)

# Test thermal response
isolates_test <- isolates_gc %>%
    filter(population == "VA") %>%
    filter(trait == "r")

mod <- lmer(value ~ temperature + site_group + temperature * site_group + (1|exp_id), data = isolates_test)
Anova(mod, type = 3)

