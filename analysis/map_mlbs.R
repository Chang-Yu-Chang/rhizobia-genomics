#' This script generate the map for MLBS and sampling site

library(tidyverse)
library(janitor)
library(cowplot)
library(ggmap)
library(ggsn) # add scale to map
source(here::here("analysis/00-metadata.R"))

dms2dec <- function(x) {
    #' convert degree-minute-second coordinate to decimal
    str_split(x, " ") %>%
        sapply(function(xx) {
            z <- as.numeric(xx[1:3])
            s <- ifelse(xx[4] %in% c("N", "E"), 1, -1)
            s*(z[1] + z[2]/60 + z[3]/3600)
        }) %>%
        return
}
site_cd <- readxl::read_excel(paste0(folder_data, "raw/High and low elevation sites_Sept 2022.xlsx"))

# Clean up the site coordinates
site_cd <- site_cd %>%
    clean_names() %>%
    mutate(latitude_dec = dms2dec(latitude_degrees_minutes_seconds)) %>%
    mutate(longitude_dec = dms2dec(longitude_degrees_minutes_seconds)) %>%
    mutate(elevation_m = elevation_feet / 3.28084) %>%
    mutate(site_group = str_sub(site,1 , 1))


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

