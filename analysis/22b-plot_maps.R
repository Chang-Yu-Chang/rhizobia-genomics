#' This script plots the map of sampling site
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

# 1. Plot both populations in one map ----
google_api <- "AIzaSyCmRu9YUWIQzq5U48yDeFPAiUYSft86pWw"
register_google(google_api)

# Philly sites
phila_center <- c(-75.25, 39.9)
phila_length <- 250
map_phila <- get_googlemap(center = phila_center, maptype = "terrain", size = c(phila_length, phila_length), style = "element:labels|visibility:off")

# MLBS site
mlbs_center <- c(-80.55, 37.35)
mlbs_length <- 250
map_mlbs <- get_googlemap(center = mlbs_center, maptype = "terrain", size = c(mlbs_length, mlbs_length), style = "element:labels|visibility:off")

# State plot
map_states <- map_data("state") %>% as_tibble %>%
    filter(region %in% c("pennsylvania", "virginia", "west virginia", "new jersey"))
state_labels <- map_states %>%
    group_by(region) %>%
    summarise(long = mean(long),lat = mean(lat)) %>%
    mutate(state = state.abb[match(region, tolower(state.name))])

p1 <- map_states %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), color = "black", fill = "NA", linewidth = .2) +
    geom_text(data = state_labels, aes(x = long, y = lat, label = state), size = 2, vjust = 1, hjust = 1) +
    # MLBS
    annotate("rect", color = rhizobia_population_colors["MLBS"], fill = NA, linewidth = 1,
             xmin = mlbs_center[1] - mlbs_length/2000, xmax = mlbs_center[1] + mlbs_length/2000,
             ymin = mlbs_center[2] - mlbs_length/2000, ymax = mlbs_center[2] + mlbs_length/2000) +
    # Phila
    annotate("rect", color = rhizobia_population_colors["Phila"], fill = NA, linewidth = 1,
             xmin = phila_center[1] - phila_length/2000, xmax = phila_center[1] + phila_length/2000,
             ymin = phila_center[2] - phila_length/2000, ymax = phila_center[2] + phila_length/2000) +
    scale_x_continuous(sec.axis = dup_axis()) +
    scale_y_continuous(sec.axis = dup_axis()) +
    coord_fixed() +
    theme_classic() +
    theme(
        plot.background =
    ) +
    labs(x = "longitude", y = "latitude")

# Phila map
p2 <- ggmap(map_phila) +
    geom_point(data = sites, aes(x = longitude_dec, y = latitude_dec, color = population, shape = population, size = population), stroke = 1, fill = NA) +
    # scalebar(x.min = map_center[1] - mlbs_length/2000, x.max = map_center[1] + mlbs_length/2000,
    #          y.min = map_center[2] - mlbs_length/2000, y.max = map_center[2] + mlbs_length/2000,
    #          dist = 5, dist_unit = "km", transform = TRUE, model = "WGS84", location = "bottomright",
    #          st.bottom = F, st.color = "black", st.size = 3, st.dist = 0.03, border.size = 0.5) +
    scale_color_manual(values = rhizobia_population_colors, name = "population") +
    scale_shape_manual(values = c(MLBS = 21, Phila = 21), name = "population") +
    scale_size_manual(values = c(MLBS = 2, Phila = 2), name = "population") +
    theme_void() +
    theme(
        plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.1, color = rhizobia_population_colors["Phila"]),
        panel.border = element_rect(color = rhizobia_population_colors["Phila"], fill = NA, linewidth = 2),
        legend.position = c(.5,.5),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(color = 1, fill = "white")
    ) +
    guides() +
    labs(x = "longitude", y = "latitude", title = "Philadelphia")

# MLBS map
p3 <- ggmap(map_mlbs) +
    geom_point(data = sites, aes(x = longitude_dec, y = latitude_dec, color = population, shape = population, size = population), stroke = 1, fill = NA) +
    # scalebar(x.min = map_center[1] - mlbs_length/2000, x.max = map_center[1] + mlbs_length/2000,
    #          y.min = map_center[2] - mlbs_length/2000, y.max = map_center[2] + mlbs_length/2000,
    #          dist = 5, dist_unit = "km", transform = TRUE, model = "WGS84", location = "bottomright",
    #          st.bottom = F, st.color = "black", st.size = 3, st.dist = 0.03, border.size = 0.5) +
    scale_color_manual(values = rhizobia_population_colors, name = "population") +
    scale_shape_manual(values = c(MLBS = 21, Phila = 21), name = "population") +
    scale_size_manual(values = c(MLBS = 2, Phila = 2), name = "population") +
    theme_void() +
    theme(
        plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.1, color = rhizobia_population_colors["MLBS"]),
        panel.border = element_rect(color = rhizobia_population_colors["MLBS"], fill = NA, linewidth = 2),
        legend.position = "none",
        legend.key.size = unit(5, "mm"), # = element_text(size = 5),
        legend.background = element_rect(color = 1, fill = "white")
    ) +
    guides() +
    labs(x = "longitude", y = "latitude", title = "MLBS")

p <- ggdraw(p1) +
    draw_plot(p2 + theme(legend.position = "none"), x = 0.6, y = 0.6, scale = .3, hjust = 0.5, vjust = 0.5) +
    draw_plot(p3, x = 0.23, y = 0.55, scale = .3, hjust = 0.5, vjust = 0.5) +
    #draw_plot(get_legend(p2), x = 0.1, y = 0.05) +
    paint_white_background()

ggsave(paste0(folder_data, "temp/22b-01-map.png"), plot = p, width = 8, height = 6)

#
sites %>%
    group_by(population) %>%
    summarize(mean_elevation = mean(elevation_m))



