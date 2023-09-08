#' This script plots the predicted t difference in the data

library(tidyverse)
library(janitor)
library(cowplot)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))


# Panel A: cartoon for hypothesis ----
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png")) + draw_text("placeholder")


# Panel B: temperature difference in sites ----
dml <- read_csv(paste0(folder_data, "temp/05-dml.csv"))
site_cd <- read_csv(paste0(folder_data, "temp/05-site_cd.csv"))
diff_tmax <- read_csv(paste0(folder_data, "temp/05-diff_tmax.csv"))
diff_tmin <- read_csv(paste0(folder_data, "temp/05-diff_tmin.csv"))
tb_season <- read_csv(paste0(folder_data, "temp/05-tb_season.csv"))
tb_month <- read_csv(paste0(folder_data, "temp/05-tb_month.csv")) %>% mutate(ymonth = factor(ymonth))

site_colors <- setNames(c(brewer.pal(6, "Blues")[3:6], brewer.pal(6, "Reds")[3:6], "gold"), site_cd$site)
site_group_colors <- setNames(c(brewer.pal(6, "Blues")[6], brewer.pal(6, "Reds")[6], "gold"), c("H", "L", "S"))
season_fills <- setNames(grey(c(0,1,0,1)), c("spring", "summer", "fall", "winter"))
month_fills <- setNames(grey(rep(c(0,1),6)), 1:12)

diff_tmax2 <- diff_tmax %>%
    group_by(yday) %>%
    summarize(mean_diff_tmax = mean(diff_tmax))
tmax_mean <- mean(diff_tmax2$mean_diff_tmax)
tmax_median <- median(diff_tmax2$mean_diff_tmax)
diff_tmin2 <- diff_tmin %>%
    group_by(yday) %>%
    summarize(mean_diff_tmin = mean(diff_tmin))
tmin_mean <- mean(diff_tmin2$mean_diff_tmin)

p2 <- left_join(diff_tmax2, diff_tmin2) %>%
    pivot_longer(cols = -yday, names_to = "mean_diff", names_prefix = "mean_diff_") %>%
    ggplot() +
    geom_boxplot(aes(x = mean_diff, y = value), width = 0.5, outlier.shape = NA) +
    geom_jitter(aes(x = mean_diff, y = value), shape = 21, width = 0.2, height = 0, size = 1, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 0.5) +
    scale_y_continuous(breaks = -2:7) +
    scale_x_discrete(breaks = c("tmax", "tmin"), labels = c(expression(t[max]), expression(t[min]))) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = grey(0.8)),
        panel.grid.minor.y = element_line(color = grey(0.95)),
        axis.title.x = element_blank()
        #axis.text = element_blank(),
        #axis.ticks.x = element_blank()
    ) +
    guides() +
    labs(y = expression(mean ~ "[" ~ t ~ "("~L~")" - t~ "("~H~")" ~ "]"))



# Panel C: cartoon for overview ----
p3 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1C.png"))

p_top <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], rel_widths = c(2,1), scale = 0.9)
p <- plot_grid(p_top, p3, nrow = 2, axis = "tblr", align = "h", labels = c("", "C"), scale = 1, rel_heights = c(1,1.5)) + paint_white_background()

ggsave(here::here("plots/Fig1.png"), p, width = 8, height = 6)








