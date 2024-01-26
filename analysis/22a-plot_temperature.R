#' This script plots the temperature data

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))

dml <- read_csv(paste0(folder_data, "temp/22-dml.csv"))
sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"))
diff_tmax <- read_csv(paste0(folder_data, "temp/22-diff_tmax.csv"))
tb_season <- read_csv(paste0(folder_data, "temp/22-tb_season.csv"))
tb_month <- read_csv(paste0(folder_data, "temp/22-tb_month.csv")) %>% mutate(ymonth = factor(ymonth))

site_colors <- setNames(c(brewer.pal(6, "Blues")[3:6], brewer.pal(6, "Reds")[3:6], "gold"), site_cd$site)
site_group_colors <- setNames(c(brewer.pal(6, "Blues")[6], brewer.pal(6, "Reds")[6], "gold"), c("H", "L", "S"))
season_fills <- setNames(grey(c(0,1,0,1)), c("spring", "summer", "fall", "winter"))
month_fills <- setNames(grey(rep(c(0,1),6)), 1:12)



# 1. Plot the raw data for the years ----
p1 <- dml %>%
    filter(site_group != "S") %>%
    ggplot() +
    geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black", linetype = 1) +
    geom_line(aes(x = yday, y = tmax_deg_c, color = population, group = site)) +
    scale_color_manual(values = rhizobia_population_colors) +
    scale_fill_manual(values = month_fills) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(-26, 40), breaks = seq(-20, 40, 10)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.1, linetype = 2),
        panel.grid.minor.y = element_line(color = "black", linewidth = 0.1, linetype = 2)
    ) +
    guides(fill = "none") +
    labs(x = "day", y = expression(t[max]))

p2 <- dml %>%
    filter(site_group != "S") %>%
    ggplot() +
    geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black", linetype = 1) +
    geom_line(aes(x = yday, y = tmin_deg_c, color = population, group = site)) +
    scale_color_manual(values = rhizobia_population_colors) +
    scale_fill_manual(values = month_fills) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(-26, 40), breaks = seq(-20, 40, 10)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.1, linetype = 2),
        panel.grid.minor.y = element_line(color = "black", linewidth = 0.1, linetype = 2)
    ) +
    guides(fill = "none") +
    labs(x = "day", y = expression(t[min]))

p <- plot_grid(p1, p2, nrow = 2, labels = c("A", "B"))
ggsave(paste0(folder_data, "temp/22a-01-daily_t.png"), p, width = 10, height = 8)

# 2. Plot the H vs L temperature ----
temp <- dml %>%
    group_by(site_group, yday) %>%
    summarize(mean_tmax_deg_c = mean(tmax_deg_c), mean_tmin_deg_c = mean(tmin_deg_c))


p1 <- temp %>%
    ggplot() +
    #geom_rect(data = temp2, aes(xmin = start, xmax = end, fill = season), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black", linetype = 1) +
    geom_line(aes(x = yday, y = mean_tmax_deg_c, color = site_group)) +
    scale_color_manual(values = site_group_colors) +
    #scale_fill_manual(values = season_fills) +
    scale_fill_manual(values = month_fills) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(-26, 33), breaks = seq(-25, 35, 5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.1, linetype = 2)
    ) +
    guides(fill = "none") +
    labs(x = "day", y = "mean max daily temperature")


p2 <- temp %>%
    ggplot() +
    #geom_rect(data = temp2, aes(xmin = start, xmax = end, fill = season), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black", linetype = 1) +
    geom_line(aes(x = yday, y = mean_tmin_deg_c, color = site_group)) +
    scale_color_manual(values = site_group_colors) +
    #scale_fill_manual(values = season_fills) +
    scale_fill_manual(values = month_fills) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(-26, 33), breaks = seq(-25, 35, 5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.1, linetype = 2)
    ) +
    guides(fill = "none") +
    labs(x = "day", y = "mean min daily temperature")

p <- plot_grid(p1, p2, nrow = 2, labels = c("A", "B"))

ggsave(paste0(folder_data, "temp/22a-02-H_vs_L.png"), p, width = 10, height = 8)

# 3. Plot the resampled temperature difference ----
## Summary stat
diff_tmax %>%
    group_by(yday) %>%
    summarize(mean_diff_tmax = mean(diff_tmax)) %>%
    summarize(mean_mean_diff_tmax = mean(mean_diff_tmax))


p <- diff_tmax %>%
    #    filter(yday == 1) %>%
    ggplot() +
    geom_histogram(aes(x = diff_tmax), color = "black", fill = "white") +
    geom_vline(xintercept = 0, color = "red", linetype = 2) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/22a-03-resample_difference_HL.png"), p, width = 4, height = 4)


# 4. For each site, use the August data  ----
p <- dml %>%
    mutate(ydate = strptime(paste("2022", yday), format="%Y %j")) %>%
    mutate(ymonth = month(ydate)) %>%
    filter(ymonth %in% 7:8) %>%
    pivot_longer(cols = ends_with("deg_c"), names_to = "t_group", values_to = "deg_c") %>%
    ggplot() +
    geom_line(aes(x = yday, y = deg_c, color = site)) +
    #scale_x_continuous(limits = c(100, 200)) +
    scale_color_manual(values = site_colors) +
    facet_wrap(~t_group) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/22a-04-august_HL.png"), p, width = 8, height = 4)


