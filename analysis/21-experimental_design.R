#' This script make figures for talks, specifically experimetnal design

library(tidyverse)
library(cowplot)
library(ggmap)
library(ggsn) # add scale to map
library(janitor)
library(waffle) #remotes::install_github("hrbrmstr/waffle")
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
    mutate(site_group = str_sub(site,1 , 1)) %>%
    filter(site != "H-1") # I did not use strain from H1


# Map sampling sites
api_secret <-  "AIzaSyByU79m1Io0KDPCw-nw1SLvGoAloaR4iSE"
register_google(key = api_secret)
mlbs_cd <- c(-80.52272538547868, 37.37521035272283)
#map_center <- c(-80.3, 37.3)
map_center <- c(-80.55, 37.35)
map_length <- 250
map <- get_googlemap(center = map_center, maptype = "terrain", size = c(map_length, map_length), style = "element:labels|visibility:off")
p1 <- ggmap(map) +
    geom_point(aes(color = "MLBS", shape = "MLBS", size = "MLBS"), x = mlbs_cd[1], y = mlbs_cd[2], stroke = 1) +
    #geom_label(aes(x = mlbs_cd[1], y = mlbs_cd[2], label = "MLBS"), vjust = 1.3, hjust = -0.3, color = "black", fill = "white") +
    geom_point(data = site_cd, aes(x = longitude_dec, y = latitude_dec, color = site_group, shape = site_group, size = site_group),
               stroke = 2, fill = "grey") +
    scalebar(x.min = map_center[1] - map_length/2000, x.max = map_center[1] + map_length/2000,
             y.min = map_center[2] - map_length/2000, y.max = map_center[2] + map_length/2000,
             dist = 5, dist_unit = "km", transform = TRUE, model = "WGS84", location = "bottomright",
             st.bottom = F, st.color = "black", st.size = 3, st.dist = 0.03, border.size = 0.5) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342", MLBS = "black"), labels = c("high", "low", "MLBS"), name = "site") +
    scale_shape_manual(values = c(H = 21, L = 21, MLBS = 4), labels = c("high", "low", "MLBS"), name = "site") +
    scale_size_manual(values = c(H = 3, L = 3, MLBS = 6), labels = c("high", "low", "MLBS"), name = "site") +
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

p <- ggdraw(p1) +
    draw_plot(p2, x = 0.15, y = 0.65, width = 0.4, height = 0.4)

ggsave(paste0(folder_data, "temp/21-01-map.png"), p, width = 5, height = 5)



# Plot counties
# usmap::us_map(regions = "counties") %>%
#     as_tibble() %>%
#     filter(full %in% c("Virginia", "West Virginia")) %>%
#     filter(!hole) %>%
#     ggplot() +
#     geom_polygon(aes(x = x, y = y, group = county), color = "black", fill = "NA", linewidth = .2)

# 2. site by elevation ----
p <- site_cd %>%
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

ggsave(paste0(folder_data, "temp/21-02-site_elevation.png"), p, width = 3, height = 3)


# 3. factorial design, plain ----
# Clean up data
treatments <- treatments %>%
    mutate(strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    mutate(strain_site = factor(strain_site, c("H", "L", "control")))

treatments %>% tabyl(strain_site, plant_site, show_missing_levels = T)
treatments %>% tabyl(strain, plant_site, show_missing_levels = T)

cc <- c(H = "#0C6291", L = "#BF4342", control = "#CBD4C2")
p1 <- treatments %>%
    group_by(strain_site, strain, plant_site) %>%
    count(.drop = F) %>%
    arrange(strain_site, plant_site) %>%
    ggplot() +
    geom_waffle(aes(values = n, fill = strain_site),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = F) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rep("#CBD4C2", 3), labels = c("high", "low", "control"), breaks = c("H", "L", "control")) +
    scale_x_continuous(breaks = seq(5, 17, 5), labels = seq(50, 170, 50)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_equal() +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = "none") +
    labs()

# factorial design, by strain strains ----
p2 <- treatments %>%
    group_by(strain_site, strain) %>%
    count(.drop = F) %>%
    arrange(strain_site, strain) %>%
    mutate(strain = factor(strain)) %>%
    ggplot() +
    geom_waffle(aes(values = n, fill = strain),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = F) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Blues")[c(3,5,7)], RColorBrewer::brewer.pal(9, "Oranges")[c(3,5,7)], "#CBD4C2")) +
    scale_x_continuous(breaks = seq(5, 17, 5), labels = seq(50, 170, 50)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_equal() +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = guide_legend(byrow = T, ncol = 3)) +
    labs()

p <- plot_grid(p1, p2, ncol = 1, axis = "tblr", align = "vh")
#ggsave(paste0(folder_data, "temp/21-05-factorial_design_rhizobia.png"), p, width = 4, height = 4)

ggsave(paste0(folder_data, "temp/21-03-factorial_design.png"), p, width = 5, height = 7)

# 4. factorial design, by maternal family ----
tt <- treatments %>%
    group_by(plant) %>%
    count(.drop = F) %>%
    arrange(desc(n)) %>%
    mutate(plant = factor(plant, plant))
p <- tt %>%
    ggplot() +
    geom_waffle(aes(values = n, fill = plant),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = F) +
    scale_fill_manual(values = rep(RColorBrewer::brewer.pal(12, "Paired"),2) %>%  setNames(tt$plant)) +
    #scale_fill_manual(values = rep(grey(c(0.1, 0.8)),20) %>%  setNames(tt$plant)) +
    scale_x_continuous(breaks = seq(5, 17, 5), labels = seq(50, 170, 50)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_equal() +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = "none") +
    labs()
ggsave(paste0(folder_data, "temp/21-04-factorial_design_plant.png"), p, width = 3.5, height = 2)




# 4. factorial design, by maternal family ----
tt <- treatments %>%
    group_by(plant) %>%
    count(.drop = F) %>%
    arrange(desc(n)) %>%
    mutate(plant = factor(plant, plant))
p <- tt %>%
    ggplot() +
    geom_waffle(aes(values = n, fill = plant),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = F) +
    scale_fill_manual(values = rep(RColorBrewer::brewer.pal(12, "Paired"),2) %>%  setNames(tt$plant)) +
    #scale_fill_manual(values = rep(grey(c(0.1, 0.8)),20) %>%  setNames(tt$plant)) +
    scale_x_continuous(breaks = seq(5, 17, 5), labels = seq(50, 170, 50)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_equal() +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = "none") +
    labs()
ggsave(paste0(folder_data, "temp/21-04-factorial_design_plant.png"), p, width = 3.5, height = 2)

# 5. factorial design, by waterblock ----
p <- treatments %>%
    group_by(waterblock, strain_site, strain) %>%
    count(.drop = F) %>%
    arrange(strain_site, strain) %>%
    mutate(strain = factor(strain)) %>%
    ggplot() +
    geom_waffle(aes(values = n, fill = strain),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = F) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Blues")[c(3,5,7)], RColorBrewer::brewer.pal(9, "Oranges")[c(3,5,7)], "#CBD4C2")) +
    scale_x_continuous(breaks = seq(5, 17, 5), labels = seq(50, 170, 50), expand = c(0,0.1)) +
    scale_y_continuous(expand = c(0,0.1)) +
    facet_wrap(~waterblock, ncol = 17) +
    coord_equal() +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = guide_legend(byrow = T, ncol = 3)) +
    labs()

ggsave(paste0(folder_data, "temp/21-05-factorial_design_waterblock.png"), p, width = 6, height = 4)

