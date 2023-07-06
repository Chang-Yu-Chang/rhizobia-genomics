#'

library(tidyverse)
library(janitor)
library(cowplot)
library(factoextra) # for plotting pca eclipse
source(here::here("analysis/00-metadata.R"))

#
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L")) %>%
    mutate(strain_site = str_sub(strain, 1, 2), strain_site_group = str_sub(strain, 1, 1))


# Clean up data
treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))

subset_ensifer <- function(tb) {
    tb %>%
        left_join(select(isolates_RDP, strain, Genus)) %>%
        drop_na()
}

gc.prm <- gc.prm %>% subset_ensifer()


# Growth traits
tt <- gc.prm %>%
    select(well, strain_site_group, all_of(c("r", "lag", "maxOD"))) %>%
    drop_na()

pcobj <- tt %>%
    select(-well, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)

p1 <- fviz_pca_ind(
    pcobj,
    label = "none",
    habillage = tt$strain_site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    scale_shape_manual(values = c(H = 16, L = 17), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = c(0.2, 0.1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()


# Extended phenotypes
tt <- treatments_M %>%
    filter(strain != "control") %>%
    select(id, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()

pcobj <- tt %>%
    select(-id, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)

p2 <- fviz_pca_ind(
    pcobj,
    label = "none",
    habillage = tt$strain_site_group,
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

p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "hv", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(here::here("plots/Fig2.png"), p, width = 8, height = 4)














