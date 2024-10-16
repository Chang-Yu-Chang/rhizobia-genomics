#' This script plots the map

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(tigris) # for getting the US census data
#library(ggmap)
#library(ggsn) # add scale to map
#library(sf) # for handling the simple features
#' On MacOS, sf depends on gdal so do this in terminal
#' brew install pkg-config
#' brew install gdal
# then in R >renv::install("sf", type = "source", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
source(here::here("metadata.R"))

sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv"))

world <- ne_countries(scale = "medium", returnclass = "sf")
sites <- data.frame(longitude = c(-80.144005, -80.109), latitude = c(26.479005, 26.83))

us_states <- states()
va_counties <- counties("VA")
va_counties$geometry

mlbs_cd <- c(-80.52272538547868, 37.37521035272283)
map_center <- c(-80.3, 37.3)
map_length <- 600

#Sys.getenv("PROJ_LIB")
us_states %>%
    filter(STUSPS %in% c("VA", "WV", "PA", "MD", "NJ")) %>%
    ggplot() +
    geom_sf() +
    geom_point(data = sites, aes(x = longitude_dec, y = latitude_dec, color = site_group, shape = site_group),
               size = 2, stroke = 1, fill = "grey") +
    scale_color_manual(values = population_colors) +
    theme_bw() +
    theme() +
    guides() +
    labs()

    #geom_point(data = sites, aes(x = longitude, y = latitude), size = 4,  shape = 23, fill = "darkred") +
    #coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE)



# Set the map center coordinates
map_center <- c(-80.3, 37.3)

# Load country boundaries from Natural Earth
countries <- ne_countries(scale = "medium", returnclass = "sf")
us_states <- states(cb = TRUE)

# Create a plot centered at the specified coordinates
ggplot(data = countries) +
    geom_sf(fill = "lightblue", color = "black") +  # Fill countries and set border color
    coord_sf(xlim = c(-81, -79), ylim = c(36.5, 37.8), expand = FALSE) +  # Adjust limits
    theme_minimal() +
    ggtitle("Map Centered at Coordinates (-80.3, 37.3)") +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title


dml <- read_csv(paste0(folder_phenotypes, "sites/dml.csv"))
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv"))
diff_vars <- read_csv(paste0(folder_phenotypes, "sites/diff_vars.csv"))
tb_season <- read_csv(paste0(folder_phenotypes, "sites/tb_season.csv"))
tb_month <- read_csv(paste0(folder_phenotypes, "sites/tb_month.csv")) %>% mutate(ymonth = factor(ymonth))

# Map sampling sites
api_secret <-  "AIzaSyByU79m1Io0KDPCw-nw1SLvGoAloaR4iSE" # You need to resigter you own API at https://mapsplatform.google.com/.
register_google(key = api_secret)
mlbs_cd <- c(-80.52272538547868, 37.37521035272283)
map_center <- c(-80.3, 37.3)
map_length <- 600
map <- get_googlemap(center = map_center, maptype = "terrain", size = c(map_length, map_length), style = "element:labels|visibility:off")
p1 <- ggmap(map) +
    geom_point(aes(color = "MLBS", shape = "MLBS"), x = mlbs_cd[1], y = mlbs_cd[2], stroke = 1, size = 3) +
    #geom_label(aes(x = mlbs_cd[1], y = mlbs_cd[2], label = "MLBS"), vjust = 1.3, hjust = -0.3, color = "black", fill = "white") +
    geom_point(data = site_cd, aes(x = longitude_dec, y = latitude_dec, color = site_group, shape = site_group),
               size = 2, stroke = 1, fill = "grey") +
    scalebar(x.min = map_center[1] - map_length/2000, x.max = map_center[1] + map_length/2000,
             y.min = map_center[2] - map_length/2000, y.max = map_center[2] + map_length/2000,
             dist = 10, dist_unit = "km", transform = TRUE, model = "WGS84", location = "bottomright",
             st.bottom = F, st.color = "black", st.size = 3, st.dist = 0.03, border.size = 0.5) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342", MLBS = "gold"), labels = c("high", "low", "MLBS"), name = "site") +
    scale_shape_manual(values = c(H = 21, L = 21, MLBS = 4), labels = c("high", "low", "MLBS"), name = "site") +
    #scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(color = 1, fill = "white")
    ) +
    guides() +
    labs(x = "longitude", y = "latitude")


# Map states
map_states <- map_data("state") %>% as_tibble %>%
    filter(region %in% c("virginia", "west virginia"))
state_labels <- map_states %>%
    group_by(region) %>%
    summarise(long = mean(long),lat = mean(lat)) %>%
    mutate(state = ifelse(region == "virginia", "VA", "WV"))

p2 <- map_states %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), color = "black", fill = "NA", linewidth = .2) +
    geom_text(data = state_labels, aes(x = long, y = lat, label = state), size = 2, vjust = 1, hjust = 1) +
    annotate("rect", color = "red", fill = NA,
             xmin = map_center[1] - map_length/2000, xmax = map_center[1] + map_length/2000,
             ymin = map_center[2] - map_length/2000, ymax = map_center[2] + map_length/2000) +
    coord_fixed() +
    theme_void() +
    theme(
        plot.background = element_rect(color = "black", fill = "white", linewidth = 1)
    ) +
    labs()

#
p <- ggdraw(p1) +
    draw_plot(p2, x = 0.6, y = 0.7, width = 0.3, height = 0.3)

ggsave(paste0(folder_data, "temp/21-01-map.png"), p, width = 5, height = 5)
