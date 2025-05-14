#'

library(tidyverse)
library(cowplot)
library(ggh4x)
library(ggrepel)
library(waffle)
library(tigris)
library(sf) # for handling the simple features
library(geodata) # for getting the elevation data
library(stars) # for converting st to sf
library(ggspatial) # for scale bar
source(here::here("metadata.R"))

# Read data ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    mutate(contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")))
sites  <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>%
    filter(population != "mid elevation") %>%
    # Use sites where the isolates were from
    filter(site %in% isolates$site) %>%
    mutate(show = T) %>%
    bind_rows(tibble(population = "PA", longitude_dec = c(-75.45, -75.05), elevation_m = 0, show = F)) %>%
    bind_rows(tibble(population = "VA", longitude_dec = c(-80.75, -80.35), elevation_m = 0, show = F))
dml <- read_csv(paste0(folder_phenotypes, "sites/dml.csv"))
tb_month <- read_csv(paste0(folder_phenotypes, "sites/tb_month.csv"))
us_elev <- elevation_30s(country = "USA", path = tempdir()) # elevation data

# Make the center. This is to use the same map range for both gradients
sites_center <- sites %>%
    drop_na() %>%
    group_by(population) %>%
    summarize(lat_mean = mean(latitude_dec), lon_mean = mean(longitude_dec), lat_max = max(latitude_dec), lat_min = min(latitude_dec), lon_max = max(longitude_dec), lon_min = min(longitude_dec)) %>%
    mutate(hjust = c(0.8, 0.8)) %>%
    mutate(width_max = max(c((lon_max-lon_min)*1.8, (lat_max-lat_min)*1.8)))

# Panel A. map ----
get_elev_sf <- function (sc, us_elev) {
    scc <- unlist(sc[,-1])
    center_coord <- scc[c("lon_mean", "lat_mean")]
    edge_length <- scc["width_max"]/2
    extent_area <- ext(center_coord[1] - edge_length, center_coord[1] + edge_length, center_coord[2] - edge_length, center_coord[2] + edge_length)
    elev1 <- crop(us_elev, extent_area)
    r_points <- as.data.frame(elev1, xy = TRUE) %>% as_tibble
    sf <- st_as_sf(r_points, coords = c("x", "y"), crs = 4326)
    tbsf <- as.data.frame(sf, xy = TRUE, cells = TRUE) %>%
        as_tibble() %>%
        mutate(
            lon = st_coordinates(geometry)[,1],
            lat = st_coordinates(geometry)[,2]
        )
    return(tbsf)
}

sftb1 <- sites_center %>%
    filter(population == "PA") %>%
    get_elev_sf(us_elev) %>%
    mutate(population = "PA")

sftb2 <- sites_center %>%
    filter(population == "VA") %>%
    get_elev_sf(us_elev) %>%
    mutate(population = "VA")

midp = 500
p1 <- bind_rows(sftb1, sftb2) %>%
    mutate(USA_elv_msk = USA_elv_msk) %>%
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = USA_elv_msk)) +
    geom_point(data = sites, aes(x = longitude_dec, y = latitude_dec, color = show), fill = "white",  size = 2, shape = 21, stroke = 1) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    scale_fill_gradient2(low = "#db7272", high = "steelblue", mid = "snow", midpoint = midp, name = "elevation (m)") +
    #annotation_scale(plot_unit = "m", width_hint = .6, location = "br", line_width = .5, text_cex = .8, height = unit(3, "mm"), pad_y = unit(1, "mm")) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-82, -71, .1)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(~population, scales = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        legend.text = element_text(size = 8) ,
        legend.title = element_text(size = 8),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(8, "mm"),
        legend.box.margin = unit(c(0,0,0,0), "mm"),
        legend.position = "right",
        strip.background = element_blank(),
        panel.spacing.x = unit(3, "mm"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = guide_colorbar(barwidth = .8, barheight = 5), color = "none") +
    labs(x = "Longitude", y = "Latitude", title = "")

# Panel B elevation ----
p2 <- sites %>%
    ggplot() +
    geom_segment(aes(x = longitude_dec, xend = longitude_dec, y = 0, yend = elevation_m), alpha = .3) +
    geom_point(aes(x = longitude_dec, y = elevation_m, color = show), fill = NA, shape = 21, stroke = 1, size = 2) +
    geom_text_repel(aes(label = site, x = longitude_dec, y = elevation_m), nudge_x = 0) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-82, -71, .1)) +
    scale_y_continuous(expand = c(.1,0)) +
    facet_wrap(~population, scales = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    guides(color = "none") +
    labs(x = "Longitude", y = "Elevation (m)")
p2
# Panel C. site temperature ----
tb_tmax <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    group_by(population) %>%
    summarize(mean_tmax = mean(tmax_deg_c))

p3 <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    ggplot() +
    geom_histogram(aes(x = tmax_deg_c), position = "identity", alpha = .5, color = "black", fill = "white", binwidth = 1) +
    geom_vline(data = tb_tmax, aes(xintercept = mean_tmax), linetype = 2) +
    geom_text(data = tb_tmax, aes(label = paste0("mean=", round(mean_tmax, 2)), x = mean_tmax), y = 100, vjust = 2, hjust = 1) +
    facet_grid2(~population, axes = "all") +
    scale_x_continuous(limits = c(5, 40)) +
    coord_flip(clip = "off") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(6, "mm"),
        panel.grid.minor.x = element_blank()
    ) +
    guides() +
    labs(x = expression("Daily maximum "(degree*C)), y = "Num. of days in Jul-Sep")

# Panel D. strain composition ----
p4 <- iso %>%
    left_join(select(isolates, population, site, genome_id)) %>%
    select(population, contig_species, site) %>%
    mutate(site = factor(site, c("fp", "src", "crp", "gp", "bg", "40th", "pms", "ppf", "L1", "L2", "H4", "H3", "H2", "L3", "L4"))) %>%
    group_by(population, site, contig_species) %>%
    count() %>%
    ggplot(aes(fill = contig_species, values = n)) +
    geom_waffle(n_rows = 1, flip = T, color = "white", size = 2, radius = unit(2, "mm")) +
    scale_fill_manual(values = species_colors) +
    scale_y_continuous(breaks = c(5, 10)) +
    facet_nested(~population+site, nest_line = element_line(colour = "black")) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.clip = "off",
        strip.background = element_blank(),
        #strip.text = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        plot.margin = unit(c(0, 12, 0, 15), unit = "mm")
    ) +
    guides(fill = guide_legend(nrow = 2)) +
    labs()

# ----
p_left <- plot_grid(
    p1, p2, p3, p4,
    ncol = 1, align = "v", axis = "rl", rel_heights = c(1.5,1,1,1),
    labels = LETTERS[1:4]
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

p <- ggdraw() +
    draw_plot(p_left, width = 1, height = 1) +
    #draw_image(here::here("plots/cartoons/Fig1.png"), scale = 1) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig1.png"), p, width = 9, height = 9)

