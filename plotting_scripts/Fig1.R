#'

library(tidyverse)
library(cowplot)
library(ggh4x)
library(ggrepel) # for annotating dots
library(waffle) # for waffle plots
library(sf) # for handling the simple features
library(geodata) # for getting the elevation data
library(stars) # for converting st to sf
library(usdata) # for temperature. data
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
# Daily max in the map region
map_range <- read_csv(paste0(folder_phenotypes, "sites/map_range.csv")) %>%
    mutate(population = rep(c("PA", "VA"), each = 31*41))
dcl <- read_csv(paste0(folder_phenotypes, "sites/dcl.csv")) %>%
    left_join(map_range)

# Panel A. map ----
p1 <- dcl %>%
    filter(yday >= 182, yday <= 273) %>%
    group_by(population, grid, longitude, latitude) %>%
    summarize(tmax = max(tmax_deg_c)) %>%
    ggplot() +
    geom_tile(aes(x = longitude, y = latitude, fill = tmax)) +
    geom_point(data = sites, aes(x = longitude_dec, y = latitude_dec, color = show), fill = "white",  size = 2, shape = 21, stroke = 1) +
    scale_fill_gradient2(low = "steelblue", high = "#db7272", mid = "snow", midpoint = 32, breaks = seq(10, 40, 2), name = "Daily maximum") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-82, -71, .1)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap2(~population, scales = "free") +
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
        panel.grid.minor = element_blank(),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.border = element_rect(color = NA, fill = NA),
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
    geom_point(aes(x = longitude_dec, y = elevation_m, color = show), fill = "white", shape = 21, stroke = 1, size = 2) +
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
        panel.grid.minor = element_blank()
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

# Chisquare
x <- iso %>%
    group_by(population, contig_species) %>%
    count() %>%
    pivot_wider(names_from = population, values_from = n, values_fill = 0)

chisq.test(x[,2:3])
