#' This script plots all plant traits measured and does LMM on all plants

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))


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

p_list <- rep(list(NA), length(trait_axis_names))
for (i in 1:length(trait_axis_names)) {
    p_list[[i]] <- treatments_M %>%
        filter(strain != "control") %>%
        drop_na(names(trait_axis_names)[i]) %>%
        plot_boxplot_pair(ytrait = names(trait_axis_names)[i], trait_axis_names[[i]])
}

p <- plot_grid(plotlist = p_list, align = "hv", axis = "tblr", ncol = 3, labels = LETTERS[1:15])

ggsave(here::here("plots/FigS10.png"), p, width = 10, height = 12)
