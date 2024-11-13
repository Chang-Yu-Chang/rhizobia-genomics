#'

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
plot_states <- function (us_states) {
    us_states %>%
        filter(STUSPS %in% c("NY", "VA", "WV", "PA", "MD", "NJ", "DE", "KY", "OH", "NC", "MI", "IN", "CT", "MA", "RI")) %>%
        ggplot() +
        #geom_sf(data = elev_sf, aes(color = USA_elv_msk)) + # elevation
        geom_sf(fill = NA, color = "grey10", linewidth = .5) +
        geom_tile(data = sites_center, aes(x = lon_mean, y = lat_mean, width = width_max, height = width_max), color = "black", fill = alpha("cornsilk", 0.5), linewidth = 0.5) +
        #scale_color_gradient(low = alpha("snow", .2), high = alpha("#d6a36e", .1), name = "elevation (m)") +
        #scale_fill_manual(values = population_colors) +
        annotation_scale(location = "bl", width_hint = .2) +
        coord_sf(expand = F, xlim = c(-87, -70), ylim = c(36, 42.5)) +
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
            # axis.text = element_blank(),
            # axis.title = element_blank(),
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
        annotation_scale(width_hint = .4, location = "br", line_width = .2, text_cex = .4, height = unit(1, "mm"), pad_y = unit(1, "mm")) +
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
            #plot.margin = unit(c(3,3,3,8), "mm"),
            plot.margin = unit(c(3,3,3,3), "mm"),
            plot.title = element_text(size = 8, hjust = .5, margin = unit(c(0,0,1,0), "mm"))
        ) +
        guides(fill = "none", color = guide_colorbar(barwidth = .5, barheight = 3)) +
        labs(x = "Longitude", y = "Latitude")
}
elev_sf1 <- get_elev_sf(us_elev, sites_center, "elevation")
p2 <- plot_elev_map(elev_sf1, "elevation", 1000)

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
        scale_color_gradient2(low = "#a642bf", high = "#0cc45f", mid = "snow", midpoint = midp, name = "time to\ncities (min)") +
        scale_fill_manual(values = population_colors) +
        annotation_scale(width_hint = .4, location = "br", line_width = .2, text_cex = .4, height = unit(1, "mm"), pad_y = unit(1, "mm")) +
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
            #plot.margin = unit(c(3,3,3,8), "mm"),
            plot.margin = unit(c(3,3,3,3), "mm"),
            plot.title = element_text(size = 8, hjust = .5, margin = unit(c(0,0,1,0), "mm"))
        ) +
        guides(fill = "none", color = guide_colorbar(barwidth = .5, barheight = 3)) +
        labs(x = "Longitude", y = "Latitude")
}
tt_sf2 <- get_tt_sf(tt, sites_center, "urbanization")
p3 <- plot_tt_map(tt_sf2, "urbanization")


# 4. Map ----
zoom_polygon1 <- polygonGrob(x = c(.385,.385,.43,.43), y = c(.5,.85,.35,.325), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
zoom_polygon2 <- polygonGrob(x = c(.7,.7,.95,.725), y = c(.6,.58,.58,.6), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
p_map <- ggdraw(p1) +
    draw_grob(zoom_polygon1) +
    draw_grob(zoom_polygon2) +
    draw_plot(p2, x = .08, y = .50, width = .35, height = .35) +
    draw_plot(p3, x = .65, y = .23, width = .35, height = .35) +
    draw_text("elevation", x = .25, y = .87, size = 10, hjust = .5) +
    draw_text("urbanization", x = .82, y = .60, size = 10, hjust = .5)

# Combine ----
p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig1.png"), scale = 1) +
    #draw_plot(p_pca, x = .3, y = .48, width = .65, height = .43) +
    draw_plot(p_map, x = .35, y = .43 , width = .55, height = .55) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig1.png"), p, width = 11, height = 7)
