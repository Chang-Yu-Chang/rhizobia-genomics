#' This script generates figure 1

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(ggmap)
#library(ggsn) # add scale to map
library(ggsci)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))

sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"), show_col_types = F)
diff_vars <- read_csv(paste0(folder_data, "temp/22-diff_vars.csv"), show_col_types = F)


# Panel A. Cartoon for characterizing divergence in two populations
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png"))

# Panel B. Plot both populations in one map
google_api <- "AIzaSyCmRu9YUWIQzq5U48yDeFPAiUYSft86pWw"
register_google(google_api)
p2 <- ggdraw()
if (FALSE) {

# Sites
mlbs_center <- c(-80.55, 37.35)
mlbs_length <- 350
get_stadiamap(center = mlbs_center, size = c(mlbs_length, mlbs_length), maptype = "stamen_terrain")
map_mlbs <- get_googlemap(center = mlbs_center, maptype = "terrain", size = c(mlbs_length, mlbs_length), style = "element:labels|visibility:off")
phila_center <- c(-75.25, 39.9)
phila_length <- 350
map_phila <- get_googlemap(center = phila_center, maptype = "terrain", size = c(phila_length, phila_length), style = "element:labels|visibility:off")

# VA map
dis <- 0.125
p2_1 <- ggmap(map_mlbs) +
    geom_point(data = filter(sites, population == "VA"), aes(x = longitude_dec, y = latitude_dec, color = site_group), shape = 21, size = 3, stroke = 1, fill = NA) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(limits = mlbs_center[1] + c(-dis, dis), breaks = ) +
    scale_y_continuous(limits = mlbs_center[2] + c(-dis, dis)) +
    coord_fixed() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.3, 0.85),
        legend.margin = margin(-5,1,1,1),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs(x = expression(longitude(~degree~W)), y = expression(latitude(~degree~N)))

## PA map
p2_2 <- ggmap(map_phila) +
    geom_point(data = filter(sites, population == "PA"), aes(x = longitude_dec, y = latitude_dec, color = site_group), shape = 21, size = 3, stroke = 1, fill = NA) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(limits = phila_center[1] + c(-dis, dis), breaks = phila_center[1] + c(0.1, 0, -0.1)) +
    scale_y_continuous(limits = phila_center[2] + c(-dis, dis), breaks = phila_center[2] + c(0.05, -0.05)) +
    coord_fixed() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.15),
        legend.key.size = unit(5, "mm"),
        legend.margin = margin(-5,1,1,1),
        legend.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs(x = expression(longitude(~degree~W)), y = expression(latitude(~degree~N)))

p2 <- plot_grid(p2_1, p2_2, ncol = 1, scale = 1, align = "vh", axis = "lrtb")
}

# Panel C. cartoon
p3 <- ggdraw()

#
p_left <- plot_grid(p1, p2, nrow = 2, labels = c("A", "B"))

p <- plot_grid(p_left, p3, nrow = 1, rel_widths = c(1, 1), labels = c("", "C")) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/Fig1.png"), p, width = 10, height = 5)



