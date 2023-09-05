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

treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))

tt_M <- treatments_M %>%
    filter(strain != "control") %>%
    select(id, strain, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()

pcobj <- tt_M %>%
    select(-id, -strain, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)
df <- as_tibble(predict(pcobj)[,1:2])
p <- df %>%
    bind_cols(select(tt_M, id, strain, strain_site_group)) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = strain_site_group, shape = strain), size = 2, stroke = 1) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    #scale_shape_manual(values = setNames(c(16,17,15,16,17,15), rhizobia_strains)) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1)
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(summary(pcobj)$importance[2,1]* 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pcobj)$importance[2,2]* 100, 1), "%)"))

ggsave(here::here("plots/FigS10.png"), p, width = 5, height = 4)

## Stats
## Explanation of PC1 + PC2
summ <- summary(pcobj)
sum(summ$importance[2,1:2]) # 0.79289

# PERMANOVA test
tt_M <- treatments_M %>%
    filter(strain != "control") %>%
    select(id, strain, strain_site_group, plant, all_of(traits), -nodule_weight_mg) %>%
    drop_na()
Y <- tt_M %>%
    select(names(trait_axis_names)[-4])

set.seed(1)
adonis2(Y ~ strain, data = tt_M, strata = tt_M$plant, permutations = 999)
# adonis2(formula = Y ~ strain, data = tt_M, permutations = 999, strata = tt_M$plant)
#           Df SumOfSqs      R2      F Pr(>F)
# strain     5   0.4599 0.05927 1.2349  0.229
# Residual  98   7.2986 0.94073
# Total    103   7.7584 1.00000











