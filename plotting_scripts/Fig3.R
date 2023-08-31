#' This script plots the plant traits

library(tidyverse)
library(janitor)
library(cowplot)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(factoextra) # for plotting pca eclipse
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)


# Clean up data
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
        facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high elevation", L="low elevation"))) +
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
Anova(mod, type = 3) # Site group does not have effect on lag time
# Response: dry_weight_mg
# Chisq Df Pr(>Chisq)
# (Intercept)       40.9288  1  1.579e-10 ***
# strain_site_group  0.8439  2     0.6558


# Panel C: PCA for extended phenotypes ----
tt_M <- treatments_M %>%
    filter(strain != "control") %>%
    select(id, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()

pcobj <- tt_M %>%
    select(-id, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)

p3 <- fviz_pca_ind(
    pcobj,
    label = "none",
    habillage = tt_M$strain_site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    scale_shape_manual(values = c(H = 16, L = 17), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = c(0.8, 0.9),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()


## Stats
# MANOVA test plus random effect
mod <- manova(cbind(dry_weight_mg) ~ strain_site_group, data = filter(treatments_M, strain != "control"))
tidy(mod)

mod <- ano(dry_weight_mg ~ strain_site_group, data = filter(treatments_M, strain != "control"))
tidy(mod)


mod <- lmer(dry_weight_mg ~ strain_site_group + (1|waterblock), data = filter(treatments_M, strain != "control"))
Anova(mod, type = 3) # rhizobia strain has effect on shoot biomass


# mod <- lmer(dry_weight_mg ~ strain_site_group + (1|plant) + (1|waterblock), data = treatments_M)
# sepl <- iris$Sepal.Length
# petl <- iris$Petal.Length
manova(cbind(Sepal.Length, Petal.Length) ~ strain_site_group + (1|plant) + (1|waterblock), data = iris)
summary(res.man)



# Panel D the experiment for plant local adaptation ----
p4 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig3A.png")) + draw_text("placeholder for\ncartoon")

# Panel E the result of local adaptation ----
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
    scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
    facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high elevation rhizobia", L="low elevation rhizobia"))) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = NA)
    ) +
    guides(color = "none", fill = guide_legend(title = "plant origin")) +
    labs(x = "", y = trait_axis_names[[1]])

## Stat



p <- plot_grid(p1, p2, p3, p4, p5, NULL, nrow = 2, axis = "tblr", align = "h", labels = LETTERS[1:6], scale = 0.95, rel_widths = c(1,1.5,1.5)) + paint_white_background()

ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 8)














