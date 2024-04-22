#' This script plots the growth trait by species

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))
sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"))
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'))

# Three strains have more than one species. Remove them from it
removed_st <- group_by(isolates_contigs, genome_id) %>% count() %>% filter(n !=1)

#
isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_longer(cols = -c(exp_id, temperature), names_to = "trait") %>%
    left_join(isolates) %>%
    filter(!(genome_id %in% removed_st$genome_id)) %>%
    left_join(select(isolates_contigs, genome_id, species)) %>%
    mutate(species = ifelse(species %in% c("medicae", "meliloti"), species, "other"))

p <- isolates_gc %>%
    filter(temperature != "40c") %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    ggplot() +
    geom_point(aes(x = temperature, y = value, color = species, group = exp_id), shape = 21, size = 2, stroke = 1) +
    geom_line(aes(x = temperature, y = value, color = species, group = exp_id), alpha = 0.5) +
    scale_color_manual(values = brewer.pal(3, "Paired")) +
    facet_grid(trait ~ population, scales = "free_y") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
    ) +
    guides() +
    labs()

ggsave(here::here("plots/FigS9.png"), p, width = 8, height = 6)


# Test species effect
isolates_test <- isolates_gc %>% filter(trait == "r")
mod <- lmer(value ~ species + population + (1|population:site) + (1|exp_id), data = isolates_test)
Anova(mod, type = 3)

# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: value
# Chisq Df Pr(>Chisq)
# (Intercept) 18.1249  1  2.069e-05 ***
# species     26.9673  3  5.981e-06 ***
# population   0.0346  1     0.8525

#GrowthRate~Species + Region + (1|Region:Site) + (1|Strain), data=alldata



















