#' This script plots the plant traits

library(tidyverse)
library(janitor)
library(cowplot)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(factoextra) # for plotting pca eclipse
library(vegan) # for permanova
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)

treatments_HL <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group != "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))

tt_HL <- treatments_HL %>%
    filter(strain != "control") %>%
    select(id, plant, plant_site_group, strain, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()

pcobj <- tt_HL %>%
    select(-id, -plant, -plant_site_group, -strain, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)
df <- as_tibble(predict(pcobj)[,1:2])
p <- df %>%
    bind_cols(select(tt_HL, id, plant, plant_site_group, strain, strain_site_group)) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = strain_site_group, shape = plant_site_group), size = 2, stroke = 1) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("H rhizobia", "L rhizobia"), breaks = c("H", "L"), name = "elevation") +
    scale_shape_manual(values = c(H = 21, L = 22), labels = c("H plant", "L plant"), breaks = c("H", "L"), name = "elevation") +
    #scale_fill_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1)
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(summary(pcobj)$importance[2,1]* 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pcobj)$importance[2,2]* 100, 1), "%)"))

ggsave(here::here("plots/FigS11.png"), p, width = 5, height = 4)

## Stats
## Explanation of PC1 + PC2
summ <- summary(pcobj)
sum(summ$importance[2,1:2]) # 0.8053

# PERMANOVA test
tt_HL <- treatments_HL %>%
    filter(strain != "control") %>%
    select(id, plant, plant_site_group, strain, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()
Y <- tt_HL %>%
    select(names(trait_axis_names)[-4])

set.seed(1)
adonis2(Y ~  plant_site_group + strain_site_group , data = tt_HL, strata = tt_HL$plant, permutations = 999)
# adonis2(formula = Y ~ plant_site_group + strain_site_group, data = tt_HL, permutations = 999, strata = tt_HL$plant)
# Df SumOfSqs      R2      F Pr(>F)
# plant_site_group   1   0.4299 0.11487 5.4770  0.817
# strain_site_group  1   0.0158 0.00421 0.2008  0.817
# Residual          42   3.2965 0.88091
# Total             44   3.7421 1.00000











