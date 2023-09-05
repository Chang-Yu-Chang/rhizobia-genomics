#' This script plots the plant traits of the control (no rhizobia inoculation)
#'

library(tidyverse)
library(janitor)
library(cowplot)
library(broom)
# library(lme4) # for linear mixed-effect models
# library(car) # companion to Applied Regression
# library(factoextra) # for plotting pca eclipse
# library(vegan) # for permanova
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
set.seed(2)

p1 <- treatments %>%
    mutate(rr = ifelse(is.na(strain), "control", "others")) %>%
    ggplot() +
    geom_boxplot(aes(x = rr, y = nodule_number), outlier.shape = NA) +
    geom_jitter(aes(x = rr, y = nodule_number), shape = 21, height = 0, width = 0.3) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "", y = "number of nodules")

p2 <- treatments %>%
    filter(is.na(strain)) %>%
    mutate(waterblock = factor(waterblock)) %>%
    replace_na(list(nodule_number = 0, dry_weight_mg = 0)) %>%
    ggplot() +
    geom_jitter(aes(x = dry_weight_mg, y = nodule_number), width = 0.15, height = 0.15, shape = 21, size = 2, stroke = 1) +
    theme_classic() +
    theme(
        legend.position = c(0.2, 0.8),
        legend.background = element_rect(color = 1, fill = NA)
    ) +
    guides() +
    labs(x = "shoot biomass (mg)", y = "number of nodules")

p <- plot_grid(p1, p2, nrow = 1, labels = c("A", "B"), scale = 0.95, align = "hv", axis = "tblr") + paint_white_background()

ggsave(here::here("plots/FigS12.png"), p, width = 6, height = 3)
