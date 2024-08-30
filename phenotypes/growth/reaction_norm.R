#' This script compute the reaction norm

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))


#
gtwl <- gtw %>%
    #filter(temperature %in% c("25c", "30c", "35c")) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    mutate(population = case_when(
        population == "VA" ~ "elevation",
        population == "PA" ~ "urbanization"
    )) %>%
    mutate(trait = case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "maxOD" ~ "yield [OD]"
    ))

# Compute the mean
gtwlm <- gtwl %>%
    group_by(population, temperature, trait, site_group) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value)))

# 1. Plot the reaction norm
p <- gtwl %>%
    ggplot() +
    geom_line(aes(x = temperature, y = value, group = well, color = site_group), alpha = 0.1) +
    # Mean value
    geom_ribbon(data = gtwlm, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = site_group, group = site_group), inherit.aes = FALSE, alpha = 0.2) +
    geom_point(data = gtwlm, aes(x = temperature, y = mean_value, color = site_group, group = site_group)) +
    geom_line(data = gtwlm, aes(x = temperature, y = mean_value, color = site_group, group = site_group)) +
    scale_color_manual(values = site_group_colors, name = "population") +
    scale_fill_manual(values = site_group_colors, name = "population") +
    facet_grid(trait~population, scales = "free_y", switch = "y") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        strip.placement = "outside",
        axis.title.y = element_blank()
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "phenotypes/growth/01-reaction_norm.png"), p, width = 6, height = 6)






