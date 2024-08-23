#' This script plots the reaction norm

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gts <- read_csv(paste0(folder_data, 'phenotypes_analysis/growth/gts.csv'))
isolates_gc <- gts %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_longer(cols = -c(exp_id, temperature), names_to = "trait") %>%
    #unite(trait, trait, temperature) %>%
    left_join(isolates) %>%
    filter(temperature != "40c") %>%
    mutate(temperature = str_remove(temperature, "c")) %>%
    mutate(population = ifelse(population == "VA", "mountain", "city")) %>%
    mutate(population = factor(population, c("mountain", "city")))

# Growth rate
plot_traits <- function (tra = "r", y_axis) {
    isolates_gc %>%
        filter(trait == tra) %>%
        ggplot() +
        geom_point(aes(x = temperature, y = value, color = site_group, group = exp_id), shape = 21, size = 2, stroke = 1) +
        geom_line(aes(x = temperature, y = value, color = site_group, group = exp_id), alpha = 0.5) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        scale_color_manual(values = site_group_colors, name = "") +
        facet_grid(. ~ population, scales = "free_y") +
        theme_classic() +
        theme(
            strip.background = element_rect(color = NA, fill = NA),
            panel.border = element_rect(color = "black", fill = NA),
        ) +
        guides() +
        labs(x = expression(temperature~(degree*C)), y = y_axis)
}
p1 <- plot_traits(tra = "r", y_axis = expression(r~(1/hr))) + theme(axis.title.x = element_blank()) + guides(color = "none")
p2 <- plot_traits(tra = "lag", y_axis = expression(t~(hr))) + theme(axis.title.x = element_blank())
p3 <- plot_traits(tra = "maxOD", y_axis = expression(x)) + guides(color = "none")

p <- plot_grid(p1, p2, p3, align = "vh", axis = "tb", ncol = 1)
ggsave(here::here("plots/FigS5.png"), p, width = 6, height = 6)


# Test thermal response
## Mountain
isolates_test <- isolates_gc %>%
    filter(population == "VA") %>%
    filter(trait == "r")

mod <- lmer(value ~ temperature + site_group + temperature * site_group + (1|exp_id), data = isolates_test)
Anova(mod, type = 3)

## City
isolates_test <- isolates_gc %>%
    filter(population == "PA") %>%
    filter(trait == "r")

mod <- lmer(value ~ temperature + site_group + temperature * site_group + (1|exp_id), data = isolates_test)
Anova(mod, type = 3)
