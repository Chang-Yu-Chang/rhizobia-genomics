#'

library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("analysis/00-metadata.R"))

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
        geom_boxplot(aes(x = strain_site_group, y = !!sym(ytrait)), fill = "white", outlier.size = -1, color = "black") +
        geom_point(aes(x = strain_site_group, y = !!sym(ytrait), group = strain, color = strain), shape = 21, size = 2, stroke = 1, fill = NA,
                   position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
        scale_color_manual(values = rep("black", 100)) +
        scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
        facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high\nelevation", L="low\nelevation"))) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_rect(color = NA, fill = NA),
            #strip.text = element_text(size = 10, color = "black"),
            strip.text = element_blank(),
            axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_blank(),
            legend.position = "none"
        ) +
        guides(color = "none") +
        labs(x = "", y = ylab)
}


p_list <- rep(list(NA), length(traits))
for (i in 1:length(traits)) p_list[[i]] <- plot_boxplot_pair(tt, names(trait_axis_names)[i], trait_axis_names[[i]])

p <- plot_grid(plotlist = p_list, ncol = 5, align = "hv", axis = "tblr")
ggsave(here::here("plots/FigS6.png"), p, width = 15, height = 9)

# STAT ----




































