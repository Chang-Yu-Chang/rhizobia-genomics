#' This script plots the plant traits

library(tidyverse)
library(janitor)
library(cowplot)
library(factoextra) # for plotting pca eclipse
source(here::here("analysis/00-metadata.R"))

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



# Panel A: plant biomass ----
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
tt <- treatments %>% drop_na(strain)
trait_axis_names <- c(
    "dry_weight_mg" = "shoot biomass (mg)",
    "nodule_number" = "number of nodules",
    "root_weight_mg" = "root biomass (mg)",
    "nodule_weight_mg" = "nodule biomass (mg)",
    "number_of_root_tips" = "number of root tips",
    "number_of_branch_points" = "number of branch points",
    "total_root_length_px" = "root length (px)",
    "branching_frequency_per_px" = "branching frequencing (1/px)",
    "network_area_px2" = "root area (px^2)",
    "average_diameter_px" = "average diameter (px)",
    "median_diameter_px" = "median diameter (px)",
    "maximum_diameter_px" = "maximum diameter (px)",
    "perimeter_px" = "perimeter (px)",
    "volume_px3" = "volume (px^3)",
    "surface_area_px2" = "surface area (px^2)"
)


plot_boxplot_pair <- function (tb, ytrait, ylab = "") {
    tb %>%
        ggplot() +
        geom_rect(data = tibble(strain_site_group = c("H", "L")), aes(fill = strain_site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_boxplot(aes_string(x = "strain_site_group", y = ytrait), fill = "white", outlier.size = -1, color = "black") +
        geom_point(aes_string(x = "strain_site_group", y = ytrait, group = "strain", color = "strain"), shape = 21, size = 2, stroke = 1, fill = NA,
                   position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
        scale_color_manual(values = rep("black", 100)) +
        scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
        facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high elevation", L="low elevation"))) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_rect(color = NA, fill = NA),
            strip.text = element_text(size = 10, color = "black"),
            #strip.text = element_blank(),
            axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_blank(),
            legend.position = "none"
        ) +
        guides(color = "none") +
        labs(x = "", y = ylab)
}

p1 <- plot_boxplot_pair(tt, names(trait_axis_names)[1], trait_axis_names[[1]])



# Panel C: PCA for extended phenotypes ----
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


p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "hv", labels = LETTERS[1:2], scale = 0.95) + paint_white_background()

ggsave(here::here("plots/Fig3.png"), p, width = 8, height = 4)














