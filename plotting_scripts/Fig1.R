#'

library(tidyverse)
library(cowplot)
library(tigris) # for getting the US state map
library(sf) # for handling the simple features
library(geodata) # for getting the elevation data
library(stars) # for converting st to sf
library(ggspatial) # for scale bar
library(geosphere) # for computing the distance between sites
library(grid) # add polygons
#' On MacOS, sf depends on gdal so do this in terminal
#' brew install pkg-config
#' brew install gdal
# then in R > renv::install("sf", type = "source", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
Sys.setenv(PROJ_LIB = "/opt/homebrew/Cellar/proj/9.5.0/share/proj") # for crs
source(here::here("metadata.R"))

# Read data
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>% filter(population != "mid elevation") %>%
    # Use sites where the isolates were from
    filter(site %in% isolates$site)
us_states <- states() # us state map
us_elev <- elevation_30s(country = "USA", path = tempdir()) # elevation data
fi <- rast(paste0(folder_data, "raw/plants/Annual_NLCD_FctImp_2022.tiff")) # imperious space data

# Make the center. This is to use the same map range for both gradients
sites_center <- sites %>%
    group_by(gradient) %>%
    summarize(lat_mean = mean(latitude_dec), lon_mean = mean(longitude_dec), lat_max = max(latitude_dec), lat_min = min(latitude_dec), lon_max = max(longitude_dec), lon_min = min(longitude_dec)) %>%
    mutate(hjust = c(0.8, 0.8)) %>%
    mutate(width_max = max(c((lon_max-lon_min)*1.8, (lat_max-lat_min)*1.8)))



# 1. Plot the states ----
plot_states <- function (us_states) {
    us_states %>%
        filter(STUSPS %in% c("NY", "VA", "WV", "PA", "MD", "NJ", "DE", "KY", "OH", "NC", "MI", "IN", "CT", "MA", "RI")) %>%
        ggplot() +
        geom_sf(fill = NA, color = "grey10", linewidth = .5) +
        geom_tile(data = sites_center, aes(x = lon_mean, y = lat_mean, width = width_max, height = width_max), color = "black", fill = alpha("#FFC857", 0.5), linewidth = 0.5) +
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
    r_points <- as.data.frame(elev1, xy = TRUE) %>% as_tibble
    sf <- st_as_sf(r_points, coords = c("x", "y"), crs = 4326)
    return(sf)
}
plot_elev_map <- function (elev_sf, gra, midp = 1000) {
    elev_sf %>%
        ggplot() +
        geom_sf(aes(color = USA_elv_msk), alpha = .9) +
        geom_point(data = filter(sites, gradient == gra), aes(x = longitude_dec, y = latitude_dec, fill = population), color = "black", size = 2, shape = 21, stroke = 1) +
        scale_color_gradient2(low = "#db7272", high = "steelblue", mid = "snow", midpoint = midp, name = "elevation (m)") +
        scale_fill_manual(values = population_colors) +
        scale_x_continuous(limits = c(sites_center$lon_mean[1] - sites_center$width_max[1]/2, sites_center$lon_mean[1] + sites_center$width_max[1]/2), expand = c(0,0)) +
        scale_y_continuous(limits = c(sites_center$lat_mean[1] - sites_center$width_max[1]/2, sites_center$lat_mean[1] + sites_center$width_max[1]/2), expand = c(0,0)) +
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
        labs(x = "Longitude", y = "Latitude", title = "Elevation")
}
elev_sf <- get_elev_sf(us_elev, sites_center, "elevation")
p2 <- plot_elev_map(elev_sf, "elevation", 1000)


# 3. Urbanization ----
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
p3 <- plot_fi_map(fi_sf, "urbanization", 50)

# 4. combine ----
zoom_polygon1 <- polygonGrob(x = c(.17,.42,.435,.43), y = c(.65,.608,.608,.65), gp = gpar(fill = "grey", alpha = 0.3, col = NA))
zoom_polygon2 <- polygonGrob(x = c(.648,.595,.85,.658), y = c(.77,.75,.75,.77), gp = gpar(fill = "grey", alpha = 0.3, col = NA))

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig1.png"), scale = 1) +
    draw_plot(p1, x = .01, y = .37, width = .95, height = .6) +
    draw_grob(zoom_polygon1) +
    draw_grob(zoom_polygon2) +
    draw_plot(p2, x = .1, y = .65, width = .4, height = .3, hjust = 0, vjust = 0) + #width = .3, height = .4, hjust = 0, halign = 0) +
    draw_plot(p3, x = .52, y = .45, width = .4, height = .3, hjust = 0, vjust = 0) + #width = .7, height = .4, hjust = 0, halign = 0) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig1.png"), p, width = 9, height = 7)
