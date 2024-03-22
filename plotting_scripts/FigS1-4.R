#' This script plots the difference in the daily tmax in VA sites

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))

dml <- read_csv(paste0(folder_data, "temp/22-dml.csv"))
sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"))
diff_vars <- read_csv(paste0(folder_data, "temp/22-diff_vars.csv"))
tb_season <- read_csv(paste0(folder_data, "temp/22-tb_season.csv"))
tb_month <- read_csv(paste0(folder_data, "temp/22-tb_month.csv")) %>% mutate(ymonth = factor(ymonth))

season_fills <- setNames(grey(c(0,1,0,1)), c("spring", "summer", "fall", "winter"))
month_fills <- setNames(grey(rep(c(0,1),6)), 1:12)


plot_composites <- function (dml_i, diff_var, diff_var_i, p1_title, p2_title) {
    # Panel A ----
    p1 <- dml_i %>%
        ggplot() +
        geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
        geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 1) +
        geom_line(aes(x = yday, y = {{diff_var}}, color = site, alpha = site)) +
        scale_fill_manual(values = month_fills) +
        scale_color_manual(values = site_colors, name = "site") +
        scale_alpha_manual(values = site_alphas, name = "site") +
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
            legend.box.margin = margin(0,0,10,0, unit = "mm"),
            legend.background = element_rect(color = 1, fill = NA)
        ) +
        guides(fill = "none", color = guide_legend(byrow = F, ncol = 2, override.aes = list(linewidth = 2), title.position = "left")) +
        labs(x = "", y = p1_title)

    # Panel B the bootstrapped temperature ----
    dvs <- diff_var_i %>%
        group_by(yday) %>%
        summarize(diff_var = mean(diff_var))
    dv_mean <- mean(dvs$diff_var)
    dv_rng <- range(dvs$diff_var); dv_rng <- c(floor(dv_rng[1]), ceiling(dv_rng[2]))
    p2 <- dvs %>%
        ggplot() +
        geom_rect(data = tb_month, aes(xmin = start, xmax = end, fill = ymonth), ymin = -Inf, ymax = Inf, alpha = 0.1) +
        geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 1) +
        geom_hline(yintercept = dv_mean, color = "maroon", linetype = 2, linewidth = 1) +
        geom_line(aes(x = yday, y = diff_var)) +
        scale_fill_manual(values = month_fills) +
        scale_x_reverse(expand = c(0,0), breaks = (tb_month$start + tb_month$end)/2, labels = month.abb) +
        scale_y_continuous(expand = c(0,0), limits = dv_rng, breaks = seq(dv_rng[1], dv_rng[2], 1), position = "right") +
        coord_flip() +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major.x = element_line(color = "black", linewidth = 0.1, linetype = 2)
        ) +
        guides(fill = "none") +
        labs(x = "", y = p2_title)


    # Padding legend for panel A ----
    p_legend <- get_legend(p1)

    # Panel C histogram
    p3 <- dvs %>%
        ggplot() +
        geom_histogram(aes(x = diff_var), color = "black", fill = "white") +
        geom_vline(xintercept = 0, color = "black", linetype = 2, linewidth = 1) +
        geom_vline(xintercept = dv_mean, color = "maroon", linetype = 2, linewidth = 1) +
        scale_x_continuous(expand = c(0,0), limits = dv_rng, breaks = seq(dv_rng[1], dv_rng[2], 1)) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        guides() +
        labs(x = p2_title)

    #
    p_top <- plot_grid(p1 + guides(color = "none", alpha = "none"), p2, nrow = 1, labels = c("A", "B"), align = "hv", axis = "lrtb", scale = 0.95)
    p_bottom <- plot_grid(p_legend, p3, labels = c("", "C"), rel_widths = c(1.05, 1), scale = 0.95)
    p <- plot_grid(p_top, p_bottom, nrow = 2, rel_heights = c(3,1)) + theme(plot.background = element_rect(fill = "white", color = NA))
    return(p)

}
compute_mean <- function (diff_var_i) {
    dvs <- diff_var_i %>% group_by(yday) %>% summarize(diff_var = mean(diff_var))
    return(list(
        mean = mean(dvs$diff_var),
        ttest = t.test(dvs$diff_var)
        ))
}


# VA tmax
dml_i <- dml %>% filter(site_group %in% c("high elevation", "low elevation"), site_group != "mid elevation") %>% mutate(site = factor(site, sites$site))
diff_var_i <- diff_vars %>% filter(population == "VA", variable == "tmax_deg_c")
p <- plot_composites(dml_i, tmax_deg_c, diff_var_i,
    expression(t[max](degree*C)),
    expression(mean ~ "[" ~ t[max] ~ "("~L~")" - t[max]~ "("~H~")" ~ "]"(degree*C)))
ggsave(here::here("plots/FigS1.png"), p, width = 8, height = 8)
compute_mean(diff_var_i)
# 3.080306
# t = 46.509, df = 364, p-value < 2.2e-16

# VA tmin
dml_i <- dml %>% filter(site_group %in% c("high elevation", "low elevation"), site_group != "mid elevation") %>% mutate(site = factor(site, sites$site))
diff_var_i <- diff_vars %>% filter(population == "VA", variable == "tmin_deg_c")
p <- plot_composites(dml_i, tmin_deg_c, diff_var_i,
    expression(t[min](degree*C)),
    expression(mean ~ "[" ~ t[min] ~ "("~L~")" - t[min]~ "("~H~")" ~ "]"(degree*C)))
ggsave(here::here("plots/FigS2.png"), p, width = 8, height = 8)
compute_mean(diff_var_i)
# 0.9417112
# t = 10.551, df = 364, p-value < 2.2e-16


# PA tmax
dml_i <- dml %>% filter(site_group %in% c("suburban", "urban"), site_group != "mid elevation") %>% mutate(site = factor(site, sites$site))
diff_var_i <- diff_vars %>% filter(population == "PA", variable == "tmax_deg_c")
p <- plot_composites(dml_i, tmax_deg_c, diff_var_i,
    expression(t[max](degree*C)),
    expression(mean ~ "[" ~ t[max] ~ "("~L~")" - t[max]~ "("~H~")" ~ "]"(degree*C)))
ggsave(here::here("plots/FigS3.png"), p, width = 8, height = 8)
compute_mean(diff_var_i)
# 0.376837
# t = 29.586, df = 364, p-value < 2.2e-16

# PA tmin
dml_i <- dml %>% filter(site_group %in% c("suburban", "urban"), site_group != "mid elevation") %>% mutate(site = factor(site, sites$site))
diff_var_i <- diff_vars %>% filter(population == "PA", variable == "tmin_deg_c")
p <- plot_composites(dml_i, tmin_deg_c, diff_var_i,
    expression(t[min](degree*C)),
    expression(mean ~ "[" ~ t[min] ~ "("~L~")" - t[min]~ "("~H~")" ~ "]"(degree*C)))
ggsave(here::here("plots/FigS4.png"), p, width = 8, height = 8)
compute_mean(diff_var_i)
# 0.2657312
# t = 14.825, df = 364, p-value < 2.2e-16





