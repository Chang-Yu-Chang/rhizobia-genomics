#' This script plots the map

library(tidyverse)
library(cowplot)
library(tigris) # for getting the US census data
library(sf) # for handling the simple features
library(geodata) # for getting the elevation data
library(stars) # for converting st to sf
library(ggspatial) # for scale bar
library(grid)
#' On MacOS, sf depends on gdal so do this in terminal
#' brew install pkg-config
#' brew install gdal
# then in R >renv::install("sf", type = "source", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
Sys.setenv(PROJ_LIB = "/opt/homebrew/Cellar/proj/9.5.0/share/proj") # for crs
source(here::here("metadata.R"))

# Get data
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>% filter(population != "mid elevation") %>%
    # Use sites where the isolates were from
    filter(site %in% isolates$site)
us_states <- states() # us state map

sites_center <- sites %>%
    group_by(gradient) %>%
    summarize(lat_mean = mean(latitude_dec), lon_mean = mean(longitude_dec), lat_max = max(latitude_dec), lat_min = min(latitude_dec), lon_max = max(longitude_dec), lon_min = min(longitude_dec)) %>%
    mutate(hjust = c(0.8, 0.8)) %>%
    mutate(width_max = max(c((lon_max-lon_min)*1.8, (lat_max-lat_min)*1.8)))

us_elev <- elevation_30s(country = "USA", path = tempdir()) # elevation data
tt <- travel_time("city", path = tempdir()) # travel time data from Nelson, Andy, Daniel J. Weiss, Jacob van Etten, Andrea Cattaneo, Theresa S. McMenomy, and Jawoo Koo. 2019. “A Suite of Global Accessibility Indicators.” Scientific Data 6 (1): 266.

# 1. Plot the states ----
get_elev_state_sf <- function (us_elev, lon_mean, lat_mean, lon_edge, lat_edge) {
    extent_area <- ext(lon_mean - lon_edge/2,
                       lon_mean + lon_edge/2,
                       lat_mean - lat_edge/2,
                       lat_mean + lat_edge/2)
    elev1 <- crop(us_elev, extent_area)
    r_points <- as.data.frame(elev1, xy = TRUE)
    elev_sf1 <- st_as_sf(r_points, coords = c("x", "y"), crs = 4326)
    as(st_geometry(elev_sf1), "Spatial")
    return(elev_sf1)
}
elev_sf <- get_elev_state_sf(us_elev, -78.5, 39, 17, 8)

p1 <- us_states %>%
    filter(STUSPS %in% c("NY", "VA", "WV", "PA", "MD", "NJ", "DE", "KY", "OH", "NC", "MI", "IN", "CT", "MA", "RI")) %>%
    ggplot() +
    geom_sf(data = elev_sf, aes(color = USA_elv_msk)) +
    geom_sf(fill = NA, color = "grey10", linewidth = .5) +
    geom_tile(data = sites_center, aes(x = lon_mean, y = lat_mean, width = width_max, height = width_max), color = "black", fill = alpha("cornsilk", 0.5), linewidth = 0.5) +
    scale_color_gradient(low = alpha("snow", .2), high = alpha("#d6a36e", .1), name = "elevation (m)") +
    scale_fill_manual(values = population_colors) +
    annotation_scale(location = "br", width_hint = .2) +
    coord_sf(expand = F, xlim = c(-87, -70), ylim = c(36, 42.5)) +
    theme_light() +
    theme(
        panel.background = element_rect(color = NA, fill = alpha("lightblue", .3), linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.background = element_rect(color = NA, fill = "white", linewidth = .5),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        axis.line = element_line(linewidth = .5),
        axis.ticks = element_line(color = "black", linewidth = .5),
        axis.text = element_text(color = "black")
        # axis.text = element_blank(),
        # axis.title = element_blank(),
    ) +
    guides(fill = "none", color = "none") +
    labs(x = "Longitude", y = "Latitude")

ggsave(paste0(folder_phenotypes, "sites/maps/01-states.png"), p1, width = 6, height = 5)

# 2. Elevation ----
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
plot_elev_map <- function (sff, gra, midp = 1000) {
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
            plot.background = element_rect(color = "black", fill = "cornsilk1", linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = .5),
            legend.text = element_text(size = 5) ,
            legend.title = element_text(size = 5),
            legend.key.width = unit(3, "mm"),
            legend.key.height = unit(5, "mm"),
            legend.box.margin = unit(c(0,0,0,0), "mm"),
            plot.margin = unit(c(3,3,3,8), "mm"),
            plot.title = element_text(size = 8, hjust = .5, margin = unit(c(0,0,1,0), "mm"))
        ) +
        guides(fill = "none") +
        labs(x = "Longitude", y = "Latitude", title = paste0(gra, " gradient"))
}
elev_sf1 <- get_elev_sf(us_elev, sites_center, "elevation")
p2 <- plot_elev_map(elev_sf1, "elevation", 1000)

ggsave(paste0(folder_phenotypes, "sites/maps/02-elevation.png"), p2, width = 6, height = 5)

# 3. Urbanization ----
get_tt_sf <- function (tt, sites_center, gra) {
    scc <- unlist(filter(sites_center, gradient == gra)[,-1])
    center_coord <- scc[c("lon_mean", "lat_mean")]
    edge_length <- scc["width_max"]/2
    extent_area <- ext(center_coord[1] - edge_length, center_coord[1] + edge_length, center_coord[2] - edge_length, center_coord[2] + edge_length)
    elev1 <- crop(tt, extent_area)
    r_points <- as.data.frame(elev1, xy = TRUE)
    elev_sf1 <- st_as_sf(r_points, coords = c("x", "y"), crs = 4326)
    as(st_geometry(elev_sf1), "Spatial")
    return(elev_sf1)
}
plot_tt_map <- function (sff, gra, midp = 55) {
    sff %>%
        ggplot() +
        geom_sf(aes(color = travel_time_to_cities_1)) +
        geom_point(data = filter(sites, gradient == gra), aes(x = longitude_dec, y = latitude_dec, fill = population), size = 2, shape = 21) +
        scale_color_gradient2(low = "#a642bf", high = "#0cc45f", mid = "snow", midpoint = midp, name = "time to cities (min)") +
        scale_fill_manual(values = population_colors) +
        annotation_scale(width_hint = .4, location = "br", line_width = .2, text_cex = .4, height = unit("1", "mm")) +
        coord_sf(clip = "off", expand = F) +
        theme_void() +
        theme(
            plot.background = element_rect(color = "black", fill = "cornsilk1", linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = .5),
            legend.text = element_text(size = 5) ,
            legend.title = element_text(size = 5),
            legend.key.width = unit(3, "mm"),
            legend.key.height = unit(5, "mm"),
            legend.box.margin = unit(c(0,0,0,0), "mm"),
            plot.margin = unit(c(3,3,3,8), "mm"),
            plot.title = element_text(size = 8, hjust = .5, margin = unit(c(0,0,1,0), "mm"))
        ) +
        guides(fill = "none") +
        labs(x = "Longitude", y = "Latitude", title = paste0(gra, " gradient"))
}
tt_sf2 <- get_tt_sf(tt, sites_center, "urbanization")
p3 <- plot_tt_map(tt_sf2, "urbanization")

ggsave(paste0(folder_phenotypes, "sites/maps/03-urbanization.png"), p3, width = 6, height = 5)

# 4. Combine ----
zoom_polygon1 <- polygonGrob(x = c(.342,.342,.415,.415), y = c(.55,.85,.345,.305), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
zoom_polygon2 <- polygonGrob(x = c(.70,.695,.73,.97), y = c(.6,.6,.53,.53), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
p <- ggdraw(p1) +
    draw_grob(zoom_polygon1) +
    draw_grob(zoom_polygon2) +
    draw_plot(p2, x = .08, y = .55, width = .3, height = .3) +
    draw_plot(p3, x = .7, y = .23, width = .3, height = .3)

ggsave(paste0(folder_phenotypes, "sites/maps/04-combined.png"), p, width = 8, height = 5)

# 5. Urbanization map by fraction impervious space ----
fi <- rast(paste0(folder_data, "raw/plants/Annual_NLCD_FctImp_2022.tiff")) # imperious space data

get_fi_sf <- function (fi) {
    project(fi, "EPSG:4326") %>%
        terra::as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        as_tibble %>%
        rename(fctimp = Annual_NLCD_FctImp_2022) %>%
        st_as_sf(coords = c("x", "y"), crs = 4326)

}
plot_fi_map <- function (fi_sf, gra, midp = 50) {
    fi_sf %>%
        ggplot() +
        geom_sf(aes(color = fctimp), alpha = .8) +
        geom_point(data = filter(sites, gradient == gra), aes(x = longitude_dec, y = latitude_dec, fill = population), color = "black", size = 2, shape = 21, stroke = 1) +
        scale_color_gradient2(low = "#0cc45f", high = "#a642bf", mid = "snow", midpoint = 50, name = "frac.\nimpervious\nsurface", breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
        scale_fill_manual(values = population_colors) +
        scale_x_continuous(limits = c(sites_center$lon_mean[2] - sites_center$width_max[2]/2, sites_center$lon_mean[2] + sites_center$width_max[2]/2), expand = c(0,0)) +
        scale_y_continuous(limits = c(sites_center$lat_mean[2] - sites_center$width_max[2]/2, sites_center$lat_mean[2] + sites_center$width_max[2]/2), expand = c(0,0)) +
        annotation_scale(width_hint = .6, location = "br", line_width = .5, text_cex = .8, height = unit(3, "mm"), pad_y = unit(1, "mm")) +
        coord_sf() +
        theme_void() +
        theme(
            plot.background = element_rect(color = "black", fill = "snow", linewidth = 1),
            legend.text = element_text(size = 8) ,
            legend.title = element_text(size = 8),
            legend.key.width = unit(5, "mm"),
            legend.key.height = unit(8, "mm"),
            legend.box.margin = unit(c(0,0,0,0), "mm"),
            plot.margin = unit(c(2,3,3,3), "mm"),
            plot.title = element_text(size = 8, hjust = .5, margin = unit(c(0,0,1,0), "mm"))
        ) +
        guides(fill = "none", color = guide_colorbar(barwidth = .8, barheight = 5)) +
        labs(x = "Longitude", y = "Latitude", title = "Urbanization")
}
fi_sf <- get_fi_sf(fi)
p <- plot_fi_map(fi_sf, "urbanization", 50)

ggsave(paste0(folder_phenotypes, "sites/maps/05-frac_impervious_surface.png"), p, width = 8, height = 5)
