#' Sampling site climates

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
tb_month <- read_csv(paste0(folder_phenotypes, "sites/tb_month.csv")) %>%  # month by day
    mutate(
        start = start - 0.5, end = end + 0.5,
        shade = rep(c("odd", "even"), 6),
        name = month.abb,
    ) %>%
    rowwise() %>%
    mutate(name_start = mean(c(start, end)))

# Panel A Every day ----
p1 <- dml %>%
    ggplot() +
    geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = shade), ymin = -Inf, ymax = Inf, alpha = .3) +
    geom_text(data = tb_month, aes(label = name, x = name_start, y = 37), size = 3) +
    geom_line(aes(x = yday, y = tmax_deg_c, group = site, color = population, linetype = "T max")) +
    geom_line(aes(x = yday, y = tmin_deg_c, group = site, color = population, linetype = "T min"), alpha = .5) +
    geom_rect(xmin = 182, xmax = 274, ymin = 0, ymax = 40, color = "black", fill = NA) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = c(VA = "steelblue", PA = "#db7272")) +
    scale_fill_manual(values = c(odd = "grey90", even = "grey60")) +
    scale_linetype_manual(values = c(`T max` = 1, `T min` = 2)) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "Days", y = expression("Temperature "(degree*C)))

# Panel B three months ----
p2 <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    ggplot() +
    geom_rect(data = filter(tb_month, ymonth %in% 7:9), aes(xmin = start, xmax = end, fill = shade), ymin = -Inf, ymax = Inf, alpha = .3) +
    geom_text(data = filter(tb_month, ymonth %in% 7:9), aes(label = name, x = name_start, y = 37)) +
    geom_line(aes(x = yday, y = tmax_deg_c, group = site, color = population, linetype = "T max")) +
    geom_line(aes(x = yday, y = tmin_deg_c, group = site, color = population, linetype = "T min"), alpha = .5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 40), expand = c(0,0)) +
    scale_color_manual(values = c(VA = "steelblue", PA = "#db7272")) +
    scale_fill_manual(values = c(odd = "grey90", even = "grey60")) +
    scale_linetype_manual(values = c(`T max` = 1, `T min` = 2)) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    guides(fill = "none", linetype = "none", color = "none") +
    labs(x = "Days", y = expression("Temperature "(degree*C)))

# Panel C ----
p3 <- dml %>%
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
    scale_x_continuous(breaks = seq(0, 40, 10), limits = c(0, 40), expand = c(0,0)) +
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
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        plot.background = element_rect(fill = NA, color = NA)
    ) +
    guides(fill = "none") +
    labs(x = expression("Daily maximum "(degree*C)), y = "Num. of days in Jul-Sep") +
    ggtitle("T max")


# Panel D ----
p4 <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    mutate(
        population = factor(population, c("VA", "PA")),
        site = factor(site, sites$site)
    ) %>%
    drop_na(site) %>%
    ggplot() +
    geom_histogram(aes(x = tmin_deg_c, fill = population), position = "identity", alpha = .6, color = "black", binwidth = 1) +
    scale_fill_manual(values = c(VA = "steelblue", PA = "#db7272")) +
    scale_y_continuous(breaks = seq(0, 150, 20)) +
    scale_x_continuous(breaks = seq(0, 40, 10), limits = c(0, 40), expand = c(0,0)) +
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
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        plot.background = element_rect(fill = NA, color = NA)
    ) +
    guides() +
    labs(x = expression("Daily maximum "(degree*C)), y = "Num. of days in Jul-Sep") +
    ggtitle("T min")



# ----
p <- plot_grid(
    p1,
    plot_grid(p2, p3, p4, nrow = 1, rel_widths = c(2,1,1), align = "h", axis = "tb", labels = LETTERS[2:4]),
    nrow = 2, labels = c("A", ""), scale = .95) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/FigS1.png"), p, width = 10, height = 8)

#ymd(paste0(2022, "-01-01")) + days(99)
