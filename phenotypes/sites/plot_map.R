#' This script plots the map

library(tidyverse)
library(cowplot)
library(tigris) # for getting the US census data
library(sf) # for handling the simple features
library(geodata) # for getting the elevation data
library(stars) # for converting st to sf
library(ggspatial) # for scale bar
library(grid)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(elevatr)
#' On MacOS, sf depends on gdal so do this in terminal
#' brew install pkg-config
#' brew install gdal
# then in R >renv::install("sf", type = "source", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
source(here::here("metadata.R"))

#
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>% filter(population != "mid elevation") %>%
    filter(population != "mid elevation")
us_states <- states()
sites_center <- sites %>%
    group_by(gradient) %>%
    summarize(lat_mean = mean(latitude_dec), lon_mean = mean(longitude_dec), lat_max = max(latitude_dec), lat_min = min(latitude_dec), lon_max = max(longitude_dec), lon_min = min(longitude_dec)) %>%
    mutate(hjust = c(0.8, 0.8)) %>%
    rowwise() %>%
    mutate(width_max = max(c((lon_max-lon_min)*1.8, (lat_max-lat_min)*1.8)))


us_elev <- elevation_30s(country = "USA", path = tempdir())


get_elev_sf <- function (us_elev, sites_center, gra) {
    scc <- unlist(filter(sites_center, gradient == gra)[,-1])
    center_coord <- scc[c("lon_mean", "lat_mean")]
    edge_length <- scc["width_max"]/2
    extent_area <- ext(center_coord[1] - edge_length, center_coord[1] + edge_length, center_coord[2] - edge_length, center_coord[2] + edge_length)
    elev1 <- crop(us_elev, extent_area)
    r_points <- as.data.frame(elev1, xy = TRUE)
    elev_sf1 <- st_as_sf(r_points, coords = c("x", "y"), crs = 4326)
    as(st_geometry(elev_sf1), "Spatial")
    return(elev_sf1)
}

# Plot the states
p1 <- us_states %>%
    filter(STUSPS %in% c("VA", "WV", "PA", "MD", "NJ", "DE")) %>%
    ggplot() +
    geom_sf(fill = "grey90", color = "grey10", linewidth = .5) +
    geom_tile(data = sites_center, aes(x = lon_mean, y = lat_mean, width = width_max, height = width_max), color = "black", fill = alpha("cornsilk", 0.5), linewidth = 0.5) +
    geom_point(data = sites, aes(x = longitude_dec, y = latitude_dec, fill = population), size = 2, shape = 21) +
    scale_fill_manual(values = population_colors) +
    annotation_scale(location = "tr") +
    coord_sf(clip = "off", expand = F, xlim = c(-87, -70), ylim = c(36, 42.5)) +
    theme_light() +
    theme(
        #panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        plot.background = element_rect(color = NA, fill = "white", linewidth = .5),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(linewidth = .5),
        axis.ticks = element_line(color = "black", linewidth = .5),
        axis.text = element_text(color = "black")
        # axis.text = element_blank(),
        # axis.ticks = element_blank(),
        # axis.title = element_blank(),
    ) +
    guides(fill = "none") +
    labs(x = "Longitude", y = "Latitude")

ggsave(paste0(folder_phenotypes, "sites/map-01-states.png"), p1, width = 6, height = 5)

# Elevation
sf1 <- get_elev_sf(us_elev, sites_center, "elevation")
plot_map <- function (sff, gra, midp = 1000) {
    sff %>%
        ggplot() +
        geom_sf(aes(color = USA_elv_msk)) +
        geom_point(data = filter(sites, gradient == gra), aes(x = longitude_dec, y = latitude_dec, fill = population), size = 2, shape = 21) +
        scale_color_gradient2(low = "maroon", high = "steelblue", mid = "snow", midpoint = midp, name = "elevation (m)") +
        scale_fill_manual(values = population_colors) +
        annotation_scale(width_hint = .4, location = "br", line_width = .2, text_cex = .4, height = unit("1", "mm")) +
        coord_sf(clip = "off", expand = F) +
        theme_void() +
        theme(
            #plot.background = element_blank(),
            plot.background = element_rect(color = "black", fill = "cornsilk1", linewidth = 1),
            #panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
            axis.ticks = element_line(color = "black", linewidth = .5),
            legend.text = element_text(size = 5) ,
            legend.title = element_text(size = 5),
            legend.key.width = unit(3, "mm"),
            legend.key.height = unit(5, "mm"),
            legend.box.margin = unit(c(0,0,0,0), "mm"),
            plot.margin = unit(c(3,3,3,3), "mm"),
            plot.title = element_text(size = 8, hjust = .5, margin = unit(c(0,0,1,0), "mm"))
        ) +
        guides(fill = "none") +
        labs(x = "Longitude", y = "Latitude", title = paste0(gra, " gradient"))
}
p2 <- plot_map(sf1, "elevation", 1000)
p2
ggsave(paste0(folder_phenotypes, "sites/map-02-elevation.png"), p2, width = 6, height = 5)

# Urbanization
sf2 <- get_elev_sf(us_elev, sites_center, "urbanization")
p3 <- plot_map(sf2, "urbanization", 30)

ggsave(paste0(folder_phenotypes, "sites/map-03-urbanization.png"), p3, width = 6, height = 5)

# Combine
zoom_polygon1 <- polygonGrob(x = c(.3,.3,.4,.4), y = c(.55,.85,.38,.27), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
zoom_polygon2 <- polygonGrob(x = c(.785,.785,.715,.715), y = c(.23,.52,.635,.59), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
p <- ggdraw(p1) +
    draw_grob(zoom_polygon1) +
    draw_grob(zoom_polygon2) +
    draw_plot(p2, x = .05, y = .55, width = .3, height = .3) +
    draw_plot(p3, x = .73, y = .23, width = .3, height = .3)

ggsave(paste0(folder_phenotypes, "sites/map-04-combined.png"), p, width = 8, height = 5)
