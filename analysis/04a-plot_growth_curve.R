#' This script is used for analyzing the growth curve
#' The experiment was done on 20230526
#' Read the README.txt in raw/growth_curve/ for details

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L"))

# Subset only Ensifer
subset_ensifer <- function(tb) {
    tb %>%
        left_join(select(isolates_RDP, strain, Genus)) %>%
        drop_na()
}

gc <- gc %>% subset_ensifer()
gc_summ <- gc_summ %>% subset_ensifer()
gc.prm <- gc.prm %>% subset_ensifer()
gc.prm.stat <- gc.prm.stat %>% subset_ensifer()


# 1. Raw data ----
p <- gc %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = strain, group = well)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04a-01-gc_raw.png"), p, width = 6, height = 5)


# 2. Curve by well ----
p <- gc %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = strain), linewidth = 1) +
    theme_light() +
    facet_grid(row ~ column) +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04a-02-gc_raw_grid.png"), p, width = 20, height = 15)

# 3. Average OD by strain ----
p <- gc_summ %>%
    mutate(strain = factor(strain, list_strains)) %>%
    left_join(rhizobia) %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs, color = rhizobia_site)) +
    geom_ribbon(aes(x = t, ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs, fill = rhizobia_site), alpha = 0.2) +
    facet_wrap(~strain, ncol = 5) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA)
    ) +
    #guides(color = "none", fill = "none") +
    labs(x = "time (hrs)", y = expression(OD[600]))
ggsave(paste0(folder_data, "temp/04a-03-gc_mean.png"), p, width = 8, height = 6)

# 4. Average OD by strain, overlayed ----
p <- gc_summ %>%
    mutate(strain = factor(strain, list_strains)) %>%
    left_join(rhizobia) %>%
    ggplot(aes(group = strain)) +
    geom_line(aes(x = t, y = mean_abs, color = rhizobia_site)) +
    scale_color_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    labs(x = "time (hrs)", y = expression(OD[600]))
ggsave(paste0(folder_data, "temp/04a-04-gc_mean_overlay.png"), p, width = 4, height = 4)

# 5. example of one strain ----
p <- gc_summ %>%
    mutate(strain = factor(strain, list_strains)) %>%
    filter(strain == "H1M1R1") %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs), color = "maroon") +
    geom_ribbon(aes(x = t, ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs), alpha = 0.2, fill = "maroon") +
    scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, 12), expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.04, 0.6), expand = c(0,0)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "time (hrs)", y = expression(OD[600]))

ggsave(paste0(folder_data, "temp/04a-05-gc_mean_example.png"), p, width = 4, height = 3)

# 6. plot all traits by strains ----
p <- gc.prm %>%
    mutate(strain = factor(strain, list_strains)) %>%
    pivot_longer(cols = c(r, t.r, lag, maxOD), names_to = "gc_trait") %>%
    ggplot() +
    geom_boxplot(aes(x = strain, y = value), outlier.size = 0) +
    geom_point(aes(x = strain, y = value), shape = 21) +
    facet_wrap(gc_trait ~., scales = "free", ncol = 1, strip.position = "right") +
    coord_flip() +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04a-06-gc_trait.png"), p, width = 6, height = 10)

# 7. plot all trait mean by site ----
p1 <- gc.prm.stat %>%
    ggplot(aes(x = site, y = lag, fill = site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = "lag time (hr)")

p2 <- gc.prm.stat %>%
    ggplot(aes(x = site, y = r, fill = site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank(),
        legend.position = "right"
    ) +
    guides() +
    labs(x = "", y = expression(growth~rate(h^-1)))

p3 <- gc.prm.stat %>%
    ggplot(aes(x = site, y = maxOD, fill = site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = expression(paste("max", "[", OD[600], "]")))

p <- plot_grid(p1, p2, p3, ncol = 1, axis = "tbrl", align = "hv")
ggsave(paste0(folder_data, "temp/04a-07-gc_trait_site.png"), p, width = 6, height = 5)

##
gc_h <- gc.prm.stat %>% filter(site == "H")
gc_l <- gc.prm.stat %>% filter(site == "L")

# lag
wilcox.test(gc_h$lag, gc_l$lag) # p=0.599
# r
wilcox.test(gc_h$r, gc_l$r) # p=0.4136
# maxOD
wilcox.test(gc_h$maxOD, gc_l$maxOD) # p=0.8518


# 8. plot all trait raw by site ----
p1 <- gc.prm %>%
    ggplot(aes(x = site, y = lag, fill = site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = "lag time (hr)")

p2 <- gc.prm %>%
    ggplot(aes(x = site, y = r, fill = site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank(),
        legend.position = "right"
    ) +
    guides() +
    labs(x = "", y = expression(growth~rate(h^-1)))

p3 <- gc.prm %>%
    ggplot(aes(x = site, y = maxOD, fill = site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    coord_flip() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = expression(paste("max", "[", OD[600], "]")))

p <- plot_grid(p1, p2, p3, ncol = 1, axis = "tbrl", align = "hv")
ggsave(paste0(folder_data, "temp/04a-08-gc_trait_site.png"), p, width = 6, height = 5)


