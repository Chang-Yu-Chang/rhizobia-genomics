#' This script analyzes the trait data

renv::load()
library(tidyverse)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(broom.mixed) # for tidy up lme4
source(here::here("metadata.R"))

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
#options(contrasts=c("contr.sum", "contr.poly"))

# Clean data
plants <- read_csv(paste0(folder_data, "phenotypes_analysis/symbiosis/plants.csv")) %>%
    left_join(isolates) %>%
    clean_names()
gts <- read_csv(paste0(folder_data, "phenotypes_analysis/growth/gts.csv")) %>%
    left_join(isolates) %>%
    clean_names() %>%
    mutate(temperature = factor(temperature, paste0(c(25, 30, 35, 40), "c"))) %>%
    arrange(temperature, exp_id) %>%
    select(genome_id, everything())

gcs <- read_csv(paste0(folder_data, "phenotypes_analysis/growth/gcs.csv")) %>%
    left_join(isolates) %>%
    clean_names() %>%
    mutate(temperature = ordered(temperature, paste0(c(25, 30, 35, 40), "c"))) %>%
    arrange(temperature, exp_id)


# Strain level difference
plants %>%
    filter(genome_id %in% c("g4", "g13")) %>%
    drop_na(shoot_biomass_mg)

gts %>%
    filter(genome_id %in% c("g4", "g13")) %>%
    drop_na(max_od)





######

gcs %>%
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    filter(population == "VA") %>%
    filter(t == 30, temperature == "40c") %>%
    filter(abs > 0.05)


gcs %>%
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    filter(population == "VA") %>%
    #filter(genome_id %in% c("g4", "g13")) %>%
    mutate(group = paste0(temperature, well)) %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = site_group, group = group)) +
    scale_color_manual(values = site_group_colors) +
    #scale_color_manual(values = c(g4="#0C6291", g13="#BF4342")) +
    facet_grid(.~temperature) +
    theme_bw() +
    theme() +
    guides() +
    labs()







