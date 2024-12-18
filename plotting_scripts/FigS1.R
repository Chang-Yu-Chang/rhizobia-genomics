#' This script plots the difference in the daily tmax in elevation sites

library(tidyverse)
library(cowplot)
source(here::here("metadata.R"))

dml <- read_csv(paste0(folder_data, "phenotypes/sites/dml.csv"))
sites <- read_csv(paste0(folder_data, "phenotypes/sites/sites.csv"))
diff_vars <- read_csv(paste0(folder_data, "phenotypes/sites/diff_vars.csv"))
tb_month <- read_csv(paste0(folder_data, "phenotypes/sites/tb_month.csv")) %>% mutate(ymonth = factor(ymonth))

site_colors <- rep(c("#0C6291", "#BF4342", "#0cc45f", "#a642bf"), each = 4) %>% setNames(sites$site[-9])
site_alphas <- rep(seq(1, 0.4, length.out = 4), 4) %>% setNames(sites$site[-9])
month_fills <- setNames(grey(rep(c(0,1),6)), 1:12)

plot_composites <- function (dml_i, diff_var, diff_var_i, p1_title, p2_title, show_legend = T) {
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
            legend.position = "right",
            legend.direction = "horizontal",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.key.height = unit("10", "mm"),
            legend.box.margin = margin(0,0,10,0, unit = "mm"),
            legend.box.background = element_blank(),
            legend.background = element_rect(color = NA, fill = NA)
        ) +
        guides(fill = "none", color = guide_legend(byrow = T, nrow = 2, override.aes = list(linewidth = 2), title.position = "left")) +
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
    p_top <- plot_grid(p1 + guides(color = "none", alpha = "none"), p2, nrow = 1, labels = c("", ""), align = "hv", axis = "lrtb", scale = 0.95)
    p_bottom <- plot_grid(NULL, p3, labels = c("", ""), rel_widths = c(1.05, 1), scale = 0.95)
    if (show_legend) p_bottom <- plot_grid(p_legend, p3, labels = c("", ""), rel_widths = c(1.05, 1), scale = 0.95)
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

get_dmli <- function (gra) dml %>% filter(gradient == gra, population != "mid elevation") %>% mutate(site = factor(site, sites$site))
get_diff_vari <- function (gra, t_deg_c) diff_vars %>% filter(gradient == gra, variable == t_deg_c)


# elevation tmax
p1 <- plot_composites(
    get_dmli("elevation"), tmax_deg_c, get_diff_vari("elevation", "tmax_deg_c"),
    expression(t[a](degree*C)),
    expression(mean ~ "[" ~ t[a] ~ "("~L~")" - t[a]~ "("~H~")" ~ "]"(degree*C))
)

# elevation tmin
p2 <- plot_composites(
    get_dmli("elevation"), tmin_deg_c, get_diff_vari("elevation", "tmin_deg_c"),
    expression(t[i](degree*C)),
    expression(mean ~ "[" ~ t[i] ~ "("~L~")" - t[i]~ "("~H~")" ~ "]"(degree*C)),
    show_legend = F
)

# urbanization tmax
p3 <- plot_composites(
    get_dmli("urbanization"), tmax_deg_c, get_diff_vari("urbanization", "tmax_deg_c"),
    expression(t[a](degree*C)),
    expression(mean ~ "[" ~ t[a] ~ "("~L~")" - t[a]~ "("~H~")" ~ "]"(degree*C))
)

# urbanization tmin
p4 <- plot_composites(
    get_dmli("urbanization"), tmin_deg_c, get_diff_vari("urbanization", "tmin_deg_c"),
    expression(t[i](degree*C)),
    expression(mean ~ "[" ~ t[i] ~ "("~L~")" - t[i]~ "("~H~")" ~ "]"(degree*C)),
    show_legend = F
)

p <- plot_grid(p1, p2, p3, p4, ncol = 2, align = "h", labels = LETTERS[1:4])

ggsave(here::here("plots/FigS1.png"), p, width = 12, height = 12)


# Stat ----
compute_mean(get_diff_vari("elevation", "tmax_deg_c"))
# 2.94137
# t = 45.844, df = 364, p-value < 2.2e-16

compute_mean(get_diff_vari("elevation", "tmin_deg_c"))
# 0.8360737
# t = 9.9063, df = 364, p-value < 2.2e-16

compute_mean(get_diff_vari("urbanization", "tmax_deg_c"))
# 0.4041381
# t = 29.515, df = 364, p-value < 2.2e-16

compute_mean(get_diff_vari("urbanization", "tmin_deg_c"))
# 0.28186
# t = 14.812, df = 364, p-value < 2.2e-16

