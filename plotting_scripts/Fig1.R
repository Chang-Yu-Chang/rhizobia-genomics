#' Map

library(tidyverse)
library(cowplot)
library(ggh4x)
library(grid)       # for plotting shades
library(lme4)       # for lmer
library(car)        # for anova
library(emmeans)    # for emmeans
library(tigris)     # for getting the US state map
library(sf)         # for handling the simple features
library(stars)      # for converting st to sf
library(usdata)     # for temperature data
library(ggspatial)  # for scale bar
library(ggrepel)    # for annotating dots
library(waffle)     # for waffle plots
source(here::here("metadata.R"))

# Read data ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    mutate(contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")))
sites  <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>%
    # Use sites where the isolates were from
    filter(site %in% isolates$site)
dml <- read_csv(paste0(folder_phenotypes, "sites/dml.csv")) # daily max t at sampling sites
tb_month <- read_csv(paste0(folder_phenotypes, "sites/tb_month.csv")) # month by day
us_states <- states() # us state map
map_range <- read_csv(paste0(folder_phenotypes, "sites/map_range.csv")) %>%
    mutate(population = rep(c("PA", "VA"), each = 31*41))
dcl <- read_csv(paste0(folder_phenotypes, "sites/dcl.csv")) %>% # Daily max in the map region
    left_join(map_range)

# Panel A. map ----
## Regional map
sites_center <- sites %>%
    drop_na(site) %>%
    group_by(population) %>%
    summarize(lat_mean = mean(latitude_dec), lon_mean = mean(longitude_dec), lat_max = max(latitude_dec), lat_min = min(latitude_dec), lon_max = max(longitude_dec), lon_min = min(longitude_dec)) %>%
    mutate(hjust = c(0.8, 0.8)) %>%
    mutate(width_max = max(c((lon_max-lon_min)*1.8, (lat_max-lat_min)*1.8)))

p1 <- us_states %>%
    filter(STUSPS %in% c("NY", "VA", "WV", "PA", "MD", "NJ", "DE", "KY", "OH", "NC", "MI", "IN", "CT", "MA", "RI")) %>%
    ggplot() +
    geom_sf(fill = NA, color = "grey10", linewidth = .5) +
    geom_tile(data = sites_center, aes(x = lon_mean, y = lat_mean, width = width_max, height = width_max), color = "black", fill = alpha("#FFC857", 0.5), linewidth = 0.5) +
    annotation_scale(location = "bl", width_hint = .2) +
    coord_sf(expand = F, xlim = c(-86.5, -70), ylim = c(35, 42.5)) +
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

## Panel A inset map
plot_temp_map <- function (dcl, sites) {

    tb <- dcl %>%
        filter(yday >= 182, yday <= 273) %>%
        group_by(longitude, latitude) %>%
        summarize(tmax = max(tmax_deg_c))

    # Create a sf point object
    sf_points <- st_as_sf(tb, coords = c("longitude", "latitude"), crs = 4326)

    sf_points %>%
        ggplot() +
        geom_sf(aes(color = tmax)) +
        geom_point(data = sites, aes(x = longitude_dec, y = latitude_dec), fill = NA, size = 1, shape = 21, stroke = .6) +
        geom_text_repel(data = sites, aes(label = site, x = longitude_dec, y = latitude_dec), size = 3) +
        scale_color_gradient2(low = "steelblue", high = "#db7272", mid = "snow", midpoint = 32, breaks = seq(26, 38, 2), limits = c(26, 38), name = expression("Daily maximum "(degree*C))) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        annotation_scale(location = "bl", width_hint = .5) +
        coord_sf(expand = F) +
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
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(0,0,0,0), "mm"),
            plot.background = element_blank(),
            plot.title = element_text(hjust = 0.1, size = 8, vjust = -10, face = "bold")
        ) +
        guides(color = guide_colorbar(barwidth = 1.2, barheight = 6)) +
        labs(x = "Longitude", y = "Latitude")
}
p1_1 <- dcl %>%
    filter(population == "VA") %>%
    filter(latitude > 37.2) %>%
    plot_temp_map(filter(sites, population == "VA")) +
    labs(title = "VA")
p1_2 <- dcl %>%
    filter(population == "PA") %>%
    filter(latitude > 39.8) %>%
    plot_temp_map(filter(sites, population == "PA")) +
    labs(title = "PA")


# Panel B. site temperature ----
tb_tmax <- dml %>%
    mutate(
        population = factor(population, c("VA", "PA")),
        site = factor(site, sites$site)
    ) %>%
    filter(site %in% isolates$site) %>%
    filter(yday >= 182 & yday <= 273) %>%
    group_by(population, site) %>%
    summarize(mean_tmax = mean(tmax_deg_c))

p2 <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    mutate(
        population = factor(population, c("VA", "PA")),
        site = factor(site, sites$site)
    ) %>%
    drop_na(site) %>%
    ggplot() +
    geom_histogram(aes(x = tmax_deg_c, fill = population), position = "identity", alpha = .5, color = "black", binwidth = 1) +
    geom_vline(data = tb_tmax, aes(xintercept = mean_tmax), linewidth = .5, linetype = 2, color = "black") +
    facet_nested(~population+site, nest_line = element_line(colour = "black")) +
    scale_fill_manual(values = c(VA = "steelblue", PA = "#db7272")) +
    scale_y_continuous(breaks = c(0, 10)) +
    scale_x_continuous(breaks = seq(0, 40, 10), limits = c(5, 40), expand = c(0,0)) +
    coord_flip(clip = "off") +
    theme_classic() +
    theme(
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        panel.spacing.x = unit(2, "mm"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = .5),
        panel.grid.minor.y = element_line(color = "grey90", linewidth = .3)
    ) +
    guides(fill = "none") +
    labs(x = expression("Daily maximum "(degree*C)), y = "Num. of days in Jul-Sep of 2022")

p2_1 <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    mutate(
        population = factor(population, c("VA", "PA")),
        site = factor(site, sites$site)
    ) %>%
    drop_na(site) %>%
    ggplot() +
    geom_histogram(aes(x = tmax_deg_c, fill = population), position = "identity", alpha = .6, color = "black", binwidth = 1) +
    scale_fill_manual(values = c(VA = "steelblue", PA = "#db7272")) +
    scale_y_continuous(breaks = seq(0, 150, 20)) +
    scale_x_continuous(breaks = seq(0, 40, 10), limits = c(5, 40), expand = c(0,0)) +
    coord_flip(clip = "off") +
    theme_classic() +
    theme(
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(2, "mm"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = .5),
        panel.grid.minor.y = element_line(color = "grey90", linewidth = .3),
        legend.position = "inside",
        legend.position.inside = c(.8, .2),
        legend.background = element_rect(color = "black", fill = "white"),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        plot.background = element_rect(fill = NA, color = NA)
    ) +
    guides() +
    labs(x = expression("Daily maximum "(degree*C)), y = "")

## Stat
xx <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    mutate(
        population = factor(population, c("VA", "PA")),
        site = factor(site, sites$site)
    ) %>%
    drop_na(site) %>%
    select(population, site, tmax_deg_c, tmin_deg_c)

mod <- lmer(tmax_deg_c ~ population + (1|site), data = xx)  # daily max
Anova(mod, type = 3)
emmeans(mod, ~ population)
mod <- lmer(tmin_deg_c ~ population + (1|site), data = xx) # daily min
Anova(mod, type = 3)
emmeans(mod, ~ population)



# Panel C. strain composition ----
p3 <- iso %>%
    left_join(select(isolates, population, site, genome_id)) %>%
    select(population, contig_species, site) %>%
    mutate(
        population = factor(population, c("VA", "PA")),
        site = factor(site, sites$site)
    ) %>%
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
        panel.spacing.x = unit(2, "mm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.clip = "off",
        strip.background = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = NA),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0, 12, 0, 15), unit = "mm")
    ) +
    guides(fill = guide_legend(ncol = 1)) +
    labs()



# ----
p_main <- plot_grid(
    p1, p2, p3,
    scale = .95,
    ncol = 1, align = "v", axis = "rl", rel_heights = c(1.5,1,.8),
    labels = LETTERS[1:3]
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

p <- ggdraw() +
    draw_plot(p_main) +
    draw_grob(polygonGrob(x = c(.29,.29,.326,.326), y = c(.705,.9,.736,.718), gp = gpar(fill = "grey", alpha = 0.3, col = NA))) +
    draw_grob(polygonGrob(x = c(.54,.54,.566,.566), y = c(.86,.864,.82,.635), gp = gpar(fill = "grey", alpha = 0.3, col = NA))) +
    draw_plot(p1_1 + guides(color = "none"), x = .1, y = .7, scale = .22, halign = 0, valign = 0) +
    draw_plot(p1_2 + guides(color = "none"), x = .53, y = .62, scale = .22, halign = 0, valign = 0) +
    draw_plot(get_legend(p1_1), x = .72, y = .65, scale = .2, halign = 0, valign = 0) +
    draw_plot(p2_1, x = .74, y = .25, width = .2, height = .228, halign = 0, valign = 0) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig1.png"), p, width = 8, height = 8)

# Chisquare
x <- iso %>%
    group_by(population, contig_species) %>%
    count() %>%
    pivot_wider(names_from = population, values_from = n, values_fill = 0)

chisq.test(x[,2:3])
