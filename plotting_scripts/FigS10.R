#' This script plots the plant traits of the control (no rhizobia inoculation)

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("metadata.R"))

plants <- read_csv(paste0(folder_data, "phenotypes_analysis/symbiosis/plants.csv"))
set.seed(2)

plants_test <- plants %>%
    mutate(inoculation = case_when(
        exp_id != "control" ~ "inoculation",
        exp_id == "control" ~ "no inoculation",
    ))

p <- plants_test %>%
    ggplot(aes(x = inoculation, y = nodule_count)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 21, height = 0, width = 0.2) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "", y = "number of nodules")

ggsave(here::here("plots/FigS10.png"), p, width = 3, height = 3)

#
mod <- lmer(nodule_count ~ inoculation + (1|waterblock), data = plants_test)
Anova(mod, type = 3)

table(plants_test$inoculation)
