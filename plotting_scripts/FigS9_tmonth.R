#' This script plots the difference in the daily tmax

library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))

dml <- read_csv(paste0(folder_data, "temp/05-dml.csv"))
site_cd <- read_csv(paste0(folder_data, "temp/05-site_cd.csv"))
diff_tmax <- read_csv(paste0(folder_data, "temp/05-diff_tmax.csv"))
diff_tmin <- read_csv(paste0(folder_data, "temp/05-diff_tmin.csv"))
tb_season <- read_csv(paste0(folder_data, "temp/05-tb_season.csv"))
tb_summer <- read_csv(paste0(folder_data, "temp/05-tb_summer.csv"))
tb_month <- read_csv(paste0(folder_data, "temp/05-tb_month.csv")) %>% mutate(ymonth = factor(ymonth))

site_colors <- setNames(c(brewer.pal(6, "Blues")[3:6], brewer.pal(6, "Reds")[3:6], "gold"), site_cd$site)
site_group_colors <- setNames(c(brewer.pal(6, "Blues")[6], brewer.pal(6, "Reds")[6], "gold"), c("H", "L", "S"))
season_fills <- setNames(grey(c(0,1,0,1)), c("spring", "summer", "fall", "winter"))
month_fills <- setNames(grey(rep(c(0,1),6)), 1:12)

# Filter to only summer and only August
dml_summer <- dml %>% filter(yday >= tb_summer$start, yday <= tb_summer$end)
diff_tmax_summer <- diff_tmax %>% filter(yday >= tb_summer$start, yday <= tb_summer$end)
diff_tmin_summer <- diff_tmin %>% filter(yday >= tb_summer$start, yday <= tb_summer$end)
tb_month_summer <- tb_month %>% filter(start >= tb_summer$start, end <= tb_summer$end)

tb_month
dml_jul <- dml %>% filter(yday >= 182, yday <= 212)
diff_tmax_jul <- diff_tmax %>% filter(yday >= 182, yday <= 212)
diff_tmin_jul <- diff_tmin %>% filter(yday >= 182, yday <= 212)
tb_month_jul <- tb_month %>% filter(start >= 182, end <= 212)


# Panel A tmax in only summer
diff_tmax_summer2 <- diff_tmax_summer %>%
    group_by(yday) %>%
    summarize(mean_diff_tmax = mean(diff_tmax))

p1 <- diff_tmax_summer2 %>%
    ggplot() +
    geom_histogram(aes(x = mean_diff_tmax), color = "black", fill = "white", binwidth = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = 2, linewidth = 1) +
    scale_x_continuous(expand = c(0,0), limits = c(-3, 8), breaks = seq(-5, 10, 5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = expression(mean ~ "[" ~ t[max] ~ "("~L~")" - t[max]~ "("~H~")" ~ "]")) +
    ggtitle(expression("Growth season" ~ t[max]))

# Panel B tmin in only summer
diff_tmin_summer2 <- diff_tmin_summer %>%
    group_by(yday) %>%
    summarize(mean_diff_tmin = mean(diff_tmin))

p2 <- diff_tmin_summer2 %>%
    ggplot() +
    geom_histogram(aes(x = mean_diff_tmin), color = "black", fill = "white", binwidth = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = 2, linewidth = 1) +
    scale_x_continuous(expand = c(0,0), limits = c(-3, 8), breaks = seq(-5, 10, 5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = expression(mean ~ "[" ~ t[min] ~ "("~L~")" - t[min]~ "("~H~")" ~ "]")) +
    ggtitle(expression("Growth season" ~ t[min]))

# Panel C tmax in August
diff_tmax_jul2 <- diff_tmax_jul %>%
    group_by(yday) %>%
    summarize(mean_diff_tmax = mean(diff_tmax))

p3 <- diff_tmax_jul2 %>%
    ggplot() +
    geom_histogram(aes(x = mean_diff_tmax), color = "black", fill = "white", binwidth = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = 2, linewidth = 1) +
    scale_x_continuous(expand = c(0,0), limits = c(-3, 8), breaks = seq(-5, 10, 5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = expression(mean ~ "[" ~ t[max] ~ "("~L~")" - t[max]~ "("~H~")" ~ "]")) +
    ggtitle(expression("July" ~ t[max]))

# Panel D tmin in August
diff_tmin_jul2 <- diff_tmin_jul %>%
    group_by(yday) %>%
    summarize(mean_diff_tmin = mean(diff_tmin))

p4 <- diff_tmin_jul2 %>%
    ggplot() +
    geom_histogram(aes(x = mean_diff_tmin), color = "black", fill = "white", binwidth = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = 2, linewidth = 1) +
    scale_x_continuous(expand = c(0,0), limits = c(-3, 8), breaks = seq(-5, 10, 5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = expression(mean ~ "[" ~ t[min] ~ "("~L~")" - t[min]~ "("~H~")" ~ "]")) +
    ggtitle(expression("July" ~ t[min]))


p <- plot_grid(p1, p2, p3, p4, nrow = 2, rel_heights = c(1,1), scale = 0.95, align = "vh", axis = "lrtb", labels = LETTERS[1:4]) + paint_white_background()
ggsave(here::here("plots/FigS9.png"), p, width = 6, height = 6)


# Stat
mean(diff_tmax_summer2$mean_diff_tmax) # 3.408841
t.test(diff_tmax_summer2$mean_diff_tmax)

mean(diff_tmin_summer2$mean_diff_tmin) # 0.790187
t.test(diff_tmin_summer2$mean_diff_tmin)

mean(diff_tmax_jul2$mean_diff_tmax) # 3.377658
t.test(diff_tmax_jul2$mean_diff_tmax)

mean(diff_tmin_jul2$mean_diff_tmin) # 0.90709
t.test(diff_tmin_jul2$mean_diff_tmin)






















