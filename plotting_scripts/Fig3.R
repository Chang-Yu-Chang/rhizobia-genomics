#' This script plots the plant traits

library(tidyverse)
library(janitor)
library(cowplot)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
#library(factoextra) # for plotting pca ellipses
library(vegan) # for permanova
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)


# Clean up data ----
treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))


treatments_HL <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group != "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))


# Panel A: cartoon for methods ----
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig3A.png")) + draw_text("placeholder for\ncartoon")


# Panel B: plant biomass ----
plot_boxplot_pair <- function (tb, ytrait, ylab = "") {
    tb %>%
        ggplot() +
        geom_rect(data = tibble(strain_site_group = c("H", "L")), aes(fill = strain_site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_boxplot(aes(x = strain_site_group, y = !!sym(ytrait)), fill = "white", outlier.size = -1, color = "black") +
        geom_point(aes(x = strain_site_group, y = !!sym(ytrait), group = strain, color = strain), shape = 21, size = 2, stroke = 1, fill = NA,
                   position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
        scale_color_manual(values = rep("black", 100)) +
        scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
        facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="H rhizobia", L="L rhizobia"))) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_rect(color = NA, fill = NA),
            strip.text = element_text(size = 10, color = "black"),
            axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_blank(),
            legend.position = "none"
        ) +
        guides(color = "none") +
        labs(x = "", y = ylab)
}

p2 <- treatments_M %>%
    filter(strain != "control") %>%
    drop_na(names(trait_axis_names)[1]) %>%
    plot_boxplot_pair(ytrait = names(trait_axis_names)[1], trait_axis_names[[1]])

## Stat
## Does the rhizobia strain have effect on shoot biomass?
mod <- lmer(dry_weight_mg ~ strain + (1|strain_site_group) + (1|plant) + (1|waterblock), data = treatments_M)
Anova(mod, type = 3) # rhizobia strain has effect on shoot biomass
# Response: dry_weight_mg
# Chisq Df Pr(>Chisq)
# (Intercept)  3.1672  1    0.07513 .
# strain      13.8485  6    0.03138 *

## Does the rhizobia sites have effect on shoot biomass?
mod <- lmer(dry_weight_mg ~ strain_site_group + (1|plant) + (1|waterblock), data = treatments_M)
Anova(mod, type = 3) # Site group does not have effect on shoot biomass
# Response: dry_weight_mg
# Chisq Df Pr(>Chisq)
# (Intercept)       40.9288  1  1.579e-10 ***
# strain_site_group  0.8439  2     0.6558


# Panel C: PCA for extended phenotypes ----
tt_M <- treatments_M %>%
    filter(strain != "control") %>%
    select(id, strain, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()

pcobj <- tt_M %>%
    select(-id, -strain, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)
df <- as_tibble(predict(pcobj)[,1:2])

p3 <- df %>%
    bind_cols(select(tt_M, strain, strain_site_group)) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = strain_site_group)) +
    geom_polygon(stat = "ellipse", aes(x = PC1, y = PC2, color = strain_site_group), fill = NA, alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("H rhizobia", "L rhizobia"), breaks = c("H", "L"), name = "elevation") +
    scale_shape_manual(values = c(H = 16, L = 17), labels = c("H rhizobia", "L rhizobia"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1),
        legend.position = "right",
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.margin = margin(0,0,0,0),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = paste0("PC1 (", round(summary(pcobj)$importance[2,1]* 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pcobj)$importance[2,2]* 100, 1), "%)"))



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
adonis2(Y ~ strain_site_group, data = tt_M, strata = tt_M$plant, permutations = 999)
# adonis2(formula = Y ~ strain_site_group, data = tt_M, permutations = 999, strata = tt_M$plant)
#                    Df SumOfSqs      R2      F Pr(>F)
# strain_site_group   1   0.0278 0.00359 0.3673  0.712
# Residual          102   7.7306 0.99641
# Total             103   7.7584 1.00000


# Panel D: the experiment for plant local adaptation ----
p4 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig3A.png")) + draw_text("placeholder for\ncartoon")

# Panel E: the result of local adaptation ----
set.seed(1)
p5 <- treatments_HL %>%
    filter(strain != "control") %>%
    drop_na(dry_weight_mg) %>%
    ggplot() +
    geom_rect(data = tibble(strain_site_group = "H"), fill = "#0C6291", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
    geom_rect(data = tibble(strain_site_group = "L"), fill = "#BF4342", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
    geom_boxplot(aes(x = strain_site_group, y = dry_weight_mg, fill = plant_site_group), alpha = 1, outlier.size = -1, color = "black") +
    geom_point(aes(x = strain_site_group, y = dry_weight_mg, group = plant_site_group, color = plant_site_group), shape = 21, size = 2, stroke = 1, fill = NA,
               position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
    scale_color_manual(values = rep("black", 100)) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("H plant", "L plant"), breaks = c("H", "L")) +
    facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="H rhizobia", L="L rhizobia"))) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.background = element_rect(fill = NA),
        legend.title = element_blank()
    ) +
    guides(color = "none") +
    labs(x = "", y = trait_axis_names[[1]])

# df_boxplot <- ggplot_build(p5)$data[[3]] # Extract the x axis
# p5 + annotate("segment", x = df_boxplot$x[1], xend = df_boxplot$x[2], y = 30, yend = 30) +
#     geom_text("", label = "*", x = mean(df_boxplot$x[1:2]), y = 35)

## Stat
## Does H rhizobia perform better with H plants than with L plants?
mod <- treatments_HL %>%
    filter(strain_site_group == "H") %>%
    lmer(dry_weight_mg ~ plant_site_group + (1|plant) + (1|waterblock), data = .)
Anova(mod, type = 3) # H rhizobia performs better when with the L plants
# Response: dry_weight_mg
#                   Chisq  Df   Pr(>Chisq)
# (Intercept)      18.9198  1  1.363e-05 ***
# plant_site_group  5.5531  1    0.01845 *


## Does L rhizobia perform better with L plants than with H plants?
mod <- treatments_HL %>%
    filter(strain_site_group == "L") %>%
    lmer(dry_weight_mg ~ plant_site_group + (1|plant) + (1|waterblock), data = .)
Anova(mod, type = 3) # No difference
# Response: dry_weight_mg
# Chisq Df Pr(>Chisq)
# (Intercept)      14.0587  1  0.0001772 ***
# plant_site_group  0.7809  1  0.3768756

## Do H plants perform better with H rhizobium than with L rhizobium?
mod <- treatments_HL %>%
    filter(plant_site_group == "H") %>%
    lmer(dry_weight_mg ~ strain_site_group + (1|strain) + (1|waterblock), data = .)
Anova(mod, type = 3)
# Response: dry_weight_mg
# Chisq Df Pr(>Chisq)
# (Intercept)       40.9996  1  1.523e-10 ***
# strain_site_group  0.1325  1     0.7159

## Do L plants perform better with L rhizobium than with H rhizobium?
mod <- treatments_HL %>%
    filter(plant_site_group == "L") %>%
    lmer(dry_weight_mg ~ strain_site_group + (1|strain) + (1|waterblock), data = .)
Anova(mod, type = 3)
# Response: dry_weight_mg
# Chisq Df Pr(>Chisq)
# (Intercept)       28.6186  1  8.813e-08 ***
# strain_site_group  0.8306  1     0.3621


# Panel F: PCA for extended phenottypes for local adpataion ----
tt_HL <- treatments_HL %>%
    filter(strain != "control") %>%
    select(id, plant, plant_site_group, strain, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()

pcobj <- tt_HL %>%
    select(-id, -plant, -plant_site_group, -strain, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)
df <- as_tibble(predict(pcobj)[,1:2])
p6 <- df %>%
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
        panel.border = element_rect(fill = NA, color = 1),
        legend.position = "right",
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.margin = margin(0,0,0,0),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(summary(pcobj)$importance[2,1]* 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pcobj)$importance[2,2]* 100, 1), "%)"))


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
# plant_site_group   1   0.4299 0.11487 5.4770  0.81
# strain_site_group  1   0.0158 0.00421 0.2008  0.81
# Residual          42   3.2965 0.88091
# Total             44   3.7421 1.00000



p <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, axis = "tblr", align = "h", labels = LETTERS[1:6], scale = 0.95, rel_widths = c(1,1.5,2)) + paint_white_background()

ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 8)














