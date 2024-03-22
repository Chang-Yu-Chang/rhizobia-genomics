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
    unite(trait, trait, temperature) %>%
    left_join(isolates)
