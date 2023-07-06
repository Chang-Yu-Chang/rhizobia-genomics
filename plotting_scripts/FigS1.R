#'

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

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
site_cd <- readxl::read_excel(paste0(folder_data, "raw/High and low elevation sites_Sept 2022.xlsx"))

# strain source
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L"))

str_sub(isolates_RDP$strain, 1, 2) %>%
    table()


# 1. map for sampling site ----
# Clean up the site coordinates
site_cd <- site_cd %>%
    clean_names() %>%
    mutate(latitude_dec = dms2dec(latitude_degrees_minutes_seconds)) %>%
    mutate(longitude_dec = dms2dec(longitude_degrees_minutes_seconds)) %>%
    mutate(elevation_m = elevation_feet / 3.28084) %>%
    mutate(site_group = str_sub(site,1 , 1) %>% factor(c("H", "S", "L"))) %>%
    filter(site != "H-1") # I did not use strain from H1


# Map sampling sites
api_secret <-  "AIzaSyByU79m1Io0KDPCw-nw1SLvGoAloaR4iSE"
register_google(key = api_secret)
mlbs_cd <- c(-80.52272538547868, 37.37521035272283)
map_center <- c(-80.55, 37.35)
map_length <- 250
map <- get_googlemap(center = map_center, maptype = "terrain", size = c(map_length, map_length), style = "element:labels|visibility:off")
p1_1 <- ggmap(map) +
    geom_point(aes(color = "MLBS", shape = "MLBS", size = "MLBS"), x = mlbs_cd[1], y = mlbs_cd[2], stroke = 1) +
    geom_point(data = site_cd, aes(x = longitude_dec, y = latitude_dec, color = site_group, shape = site_group, size = site_group),
               stroke = 2, fill = "grey") +
    scalebar(x.min = map_center[1] - map_length/2000, x.max = map_center[1] + map_length/2000,
             y.min = map_center[2] - map_length/2000, y.max = map_center[2] + map_length/2000,
             dist = 5, dist_unit = "km", transform = TRUE, model = "WGS84", location = "bottomright",
             st.bottom = F, st.color = "black", st.size = 3, st.dist = 0.03, border.size = 0.5) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342", S = "grey30", MLBS = "black"), labels = c("high", "mid", "low", "MLBS"), breaks = c("H", "S", "L", "MLBS"), name = "site") +
    scale_shape_manual(values = c(H = 21, L = 21, S = 21, MLBS = 4), labels = c("high", "mid", "low", "MLBS"), breaks = c("H", "S", "L", "MLBS"), name = "site") +
    scale_size_manual(values = c(H = 3, L = 3, S = 3, MLBS = 5), labels = c("high", "mid", "low", "MLBS"), breaks = c("H", "S", "L", "MLBS"), name = "site") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = c(0.85, 0.80),
        legend.key.size = unit(5, "mm"), # = element_text(size = 5),
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

p1_2 <- map_states %>%
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

p1 <- ggdraw(p1_1) +
    draw_plot(p1_2, x = 0.15, y = 0.65, width = 0.4, height = 0.4)


p2 <- site_cd %>%
    ggplot() +
    geom_col(aes(x = site, y = elevation_m, color = site_group), fill = "grey", linewidth = 1, width = .8) +
    #geom_col(aes(x = site, y = elevation_m), color = "black", fill = "black") +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342", MLBS = "gold"), labels = c("high", "low", "MLBS"), name = "site") +
    # scale_fill_manual(values = c(H = "#0C6291", L = "#BF4342", MLBS = "gold"), labels = c("high", "low", "MLBS"), name = "site") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1200), breaks = seq(0, 1200, 200)) +
    theme_classic() +
    theme() +
    guides(color = "none") +
    labs(x = "site", y = "elevation (m)")

p <- plot_grid(p1, p2, nrow = 1, labels = c("A", "B"), scale = c(1, 0.9)) + paint_white_background()

ggsave(here::here("plots/FigS1.png"), p, width = 8, height = 4)

