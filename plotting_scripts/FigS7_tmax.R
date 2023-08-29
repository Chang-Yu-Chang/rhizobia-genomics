#' This script plots the difference in the daily tmax

library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))

dml <- read_csv(paste0(folder_data, "temp/05-dml.csv"))
site_cd <- read_csv(paste0(folder_data, "temp/05-site_cd.csv"))
diff_tmax <- read_csv(paste0(folder_data, "temp/05-diff_tmax.csv"))
tb_season <- read_csv(paste0(folder_data, "temp/05-tb_season.csv"))
tb_month <- read_csv(paste0(folder_data, "temp/05-tb_month.csv")) %>% mutate(ymonth = factor(ymonth))

site_colors <- setNames(c(brewer.pal(6, "Blues")[3:6], brewer.pal(6, "Reds")[3:6], "gold"), site_cd$site)
site_group_colors <- setNames(c(brewer.pal(6, "Blues")[6], brewer.pal(6, "Reds")[6], "gold"), c("H", "L", "S"))
season_fills <- setNames(grey(c(0,1,0,1)), c("spring", "summer", "fall", "winter"))
month_fills <- setNames(grey(rep(c(0,1),6)), 1:12)


# Panel A the maximum temperature ----
p1 <- dml %>%
    filter(site_group != "S") %>%
    ggplot() +
    geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 1) +
    geom_line(aes(x = yday, y = tmax_deg_c, color = site)) +
    scale_color_manual(values = site_colors) +
    scale_fill_manual(values = month_fills) +
    scale_x_reverse(expand = c(0,0), breaks = (tb_month$start + tb_month$end)/2, labels = month.abb) +
    scale_y_continuous(expand = c(0,0), limits = c(-20, 35), breaks = seq(-20, 30, 10), position = "right") +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(color = "black", linewidth = 0.1, linetype = 2),
        panel.grid.minor.x = element_line(color = "black", linewidth = 0.1, linetype = 2),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.height = unit("10", "mm"),
        legend.box.margin = margin(0,0,10,0, unit = "mm")
    ) +
    guides(fill = "none", color = guide_legend(byrow = F, ncol = 2, override.aes = list(linewidth = 2))) +
    labs(x = "", y = expression(t[max]))

# Panel B the resampled temperature ----
diff_tmax2 <- diff_tmax %>%
    group_by(yday) %>%
    summarize(mean_diff_tmax = mean(diff_tmax))
p2 <- diff_tmax2 %>%
    ggplot() +
    geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 1) +
    geom_line(aes(x = yday, y = mean_diff_tmax)) +
    scale_color_manual(values = site_colors) +
    scale_fill_manual(values = month_fills) +
    scale_x_reverse(expand = c(0,0), breaks = (tb_month$start + tb_month$end)/2, labels = month.abb) +
    scale_y_continuous(expand = c(0,0), limits = c(-3, 8), breaks = seq(-5, 10, 1), position = "right") +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(color = "black", linewidth = 0.1, linetype = 2)
    ) +
    guides(fill = "none") +
    labs(x = "", y = expression(mean ~ "[" ~ t[max] ~ "("~L~")" - t[max]~ "("~H~")" ~ "]"))


# Padding legend for panel A ----
p_legend <- get_legend(p1)

# Panel C histogram
p3 <- diff_tmax2 %>%
    ggplot() +
    geom_histogram(aes(x = mean_diff_tmax), color = "black", fill = "white") +
    geom_vline(xintercept = 0, color = "black", linetype = 2, linewidth = 1) +
    scale_x_continuous(expand = c(0,0), limits = c(-3, 8), breaks = seq(-5, 10, 5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = expression(mean ~ "[" ~ t[max] ~ "("~L~")" - t[max]~ "("~H~")" ~ "]"))

p_top <- plot_grid(p1 + guides(color = "none"), p2, nrow = 1, labels = c("A", "B"), align = "hv", axis = "lrtb", scale = 0.95)
p_bottom <- plot_grid(p_legend, p3, labels = c("", "C"), rel_widths = c(1.05, 1), scale = 0.95)
p <- plot_grid(p_top, p_bottom, nrow = 2, rel_heights = c(3,1)) + paint_white_background()
ggsave(here::here("plots/FigS7.png"), p, width = 8, height = 8)


# Stat
mean(diff_tmax2$mean_diff_tmax) # 3.098698
t.test(diff_tmax2$mean_diff_tmax)






