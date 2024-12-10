#'

library(tidyverse)
library(cowplot)
library(tigris) # for getting the US census data and state map
library(sf) # for handling the simple features
library(geodata) # for getting the elevation data
library(stars) # for converting st to sf
library(ggspatial) # for scale bar
library(geosphere) # for computing the distance between sites
library(ggspatial)
library(grid)
#' On MacOS, sf depends on gdal so do this in terminal
#' brew install pkg-config
#' brew install gdal
# then in R > renv::install("sf", type = "source", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
Sys.setenv(PROJ_LIB = "/opt/homebrew/Cellar/proj/9.5.0/share/proj")
source(here::here("metadata.R"))

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
#tt <- travel_time("city", path = tempdir()) # travel time data from Nelson, Andy, Daniel J. Weiss, Jacob van Etten, Andrea Cattaneo, Theresa S. McMenomy, and Jawoo Koo. 2019. “A Suite of Global Accessibility Indicators.” Scientific Data 6 (1): 266.
city_hall <- c(-75.1652, 39.9526) # longitude, latitude
grid_size <- 0.01 # approximately 1 km
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
plot_states <- function (us_states) {
    us_states %>%
        filter(STUSPS %in% c("NY", "VA", "WV", "PA", "MD", "NJ", "DE", "KY", "OH", "NC", "MI", "IN", "CT", "MA", "RI")) %>%
        ggplot() +
        geom_sf(fill = NA, color = "grey10", linewidth = .5) +
        geom_tile(data = sites_center, aes(x = lon_mean, y = lat_mean, width = width_max, height = width_max), color = "black", fill = alpha("cornsilk", 0.5), linewidth = 0.5) +
        annotation_scale(location = "bl", width_hint = .2, ) +
        coord_sf(expand = F, xlim = c(-87, -70), ylim = c(35, 42.5)) +
        theme_light() +
        theme(
            panel.background = element_rect(color = NA, fill = NA, linewidth = 1),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.background = element_rect(color = NA, fill = NA, linewidth = .5),
            legend.position = "bottom",
            panel.grid.major = element_blank(),
            axis.line = element_line(linewidth = .5),
            axis.ticks = element_line(color = "black", linewidth = .5),
            axis.text = element_text(color = "black")
        ) +
        guides(fill = "none", color = "none") +
        labs(x = "Longitude", y = "Latitude")
}
p1 <- plot_states(us_states)

# 2. Elevation gradient ----
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
        annotation_scale(width_hint = .6, location = "br", line_width = .5, text_cex = .8, height = unit(3, "mm"), pad_y = unit(1, "mm")) +
        coord_sf(clip = "off", expand = F) +
        theme_void() +
        theme(
            plot.background = element_rect(color = "black", fill = "snow", linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = .5),
            legend.text = element_text(size = 8) ,
            legend.title = element_text(size = 8),
            legend.key.width = unit(5, "mm"),
            legend.key.height = unit(8, "mm"),
            legend.box.margin = unit(c(0,0,0,0), "mm"),
            plot.margin = unit(c(2,3,3,3), "mm"),
            plot.title = element_text(size = 8, hjust = .5, margin = unit(c(0,0,1,0), "mm"))
        ) +
        guides(fill = "none", color = guide_colorbar(barwidth = .8, barheight = 5)) +
        labs(x = "Longitude", y = "Latitude", title = "Elevation")
}
elev_sf1 <- get_elev_sf(us_elev, sites_center, "elevation")
p2 <- plot_elev_map(elev_sf1, "elevation", 1000)

# 3. Urbanization ----
get_tt_sf <- function (sites_center, grid_size = 0.001) {
    city_hall <- c(-75.1652, 39.9526) # longitude, latitude
    #grid_size <- 0.001 # approximately 100m

    dist_to_ch <- function (grid_points) {
        xx <- apply(grid_points, 1, function(row) {
            distHaversine(c(row['longitude'], row['latitude']), city_hall)
        })

        return(mutate(grid_points, distance_to_city_hall = xx))
    }

    sites_center <- sites_center %>%
        mutate(
            lat_range = map2(lat_mean-width_max/2, lat_mean+width_max/2, ~seq(.x, .y, by = grid_size)),
            lon_range = map2(lon_mean-width_max/2, lon_mean+width_max/2, ~seq(.x, .y, by = grid_size)),
            grid_points = map2(lon_range, lat_range, ~expand.grid(longitude = .x, latitude = .y)),
            gps = map(grid_points, dist_to_ch)
        )

    sc <- sites_center %>%
        select(gradient, gps) %>%
        unnest(gps) %>%
        mutate(distance_to_city_hall_km = distance_to_city_hall / 1000)
    tt_sf <- st_as_sf(sc, coords = c("longitude", "latitude"), crs = 4326)

    return(tt_sf)
}
plot_tt_map <- function (tt_sf, gra, midp = 10) {
    #gra = "urbanization"
    tt_sf %>%
        filter(gradient == gra) %>%
        ggplot() +
        geom_sf(aes(color = distance_to_city_hall_km)) +
        geom_point(data = filter(sites, gradient == gra), aes(x = longitude_dec, y = latitude_dec, fill = population), size = 2, shape = 21) +
        scale_color_gradient2(low = "#a642bf", high = "#0cc45f", mid = "snow", midpoint = midp, name = "distance to\ncity hall (km)") +
        scale_fill_manual(values = population_colors) +
        annotation_scale(width_hint = .6, location = "br", line_width = .5, text_cex = .8, height = unit(3, "mm"), pad_y = unit(1, "mm")) +
        coord_sf(clip = "off", expand = F) +
        theme_void() +
        theme(
            plot.background = element_rect(color = "black", fill = "snow", linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = .5),
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
tt_sf <- get_tt_sf(sites_center)
p3 <- plot_tt_map(tt_sf, "urbanization")


# 4. combine ----
# zoom_polygon1 <- polygonGrob(x = c(.385,.385,.43,.43), y = c(.5,.85,.35,.325), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
# zoom_polygon2 <- polygonGrob(x = c(.7,.7,.95,.725), y = c(.6,.58,.58,.6), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig1.png"), scale = 1) +
    draw_plot(p1, x = .01, y = .37, width = .95, height = .6) +
    #draw_plot(p_map, x = .35, y = .43, width = .55, height = .55) +
    # draw_grob(zoom_polygon1) +
    # draw_grob(zoom_polygon2) +
    draw_plot(p2, x = .1, y = .65, width = .4, height = .3, hjust = 0, vjust = 0) + #width = .3, height = .4, hjust = 0, halign = 0) +
    draw_plot(p3, x = .52, y = .45, width = .4, height = .3, hjust = 0, vjust = 0) + #width = .7, height = .4, hjust = 0, halign = 0) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig1.png"), p, width = 9, height = 7)
