#' This script plots the root traits

library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

# growth curves
gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L"))

# common garden
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
treatments_long <- read_csv(paste0(folder_data, "temp/11-treatments_long.csv"), show_col_types = F)
treatments_scaled <- read_csv(paste0(folder_data, "temp/11-treatments_scaled.csv"), show_col_types = F)
treatments_scaled_long <- read_csv(paste0(folder_data, "temp/11-treatments_scaled_long.csv"), show_col_types = F)

# 0. combine the data from treatments and gc ----
# Subset only Ensifer -----
subset_ensifer <- function(tb) {
    tb %>%
        left_join(select(isolates_RDP, strain, Genus)) %>%
        drop_na()
}

gc <- gc %>% subset_ensifer()
gc_summ <- gc_summ %>% subset_ensifer()
gc.prm <- gc.prm %>% subset_ensifer()
gc.prm.stat <- gc.prm.stat %>% subset_ensifer()

# Combine the data from rhizobia and ----
treatments_long_stat <- treatments_long %>%
    filter(trait %in% c("dry_weight_mg", "nodule_number", "root_weight_mg", "total_root_length_px", "branching_frequency_per_px", "network_area_px2", "average_diameter_px")) %>%
    group_by(rhizobia_site, rhizobia, trait) %>%
    drop_na(rhizobia, value) %>%
    summarize(trait_mean = mean(value, na.rm = T), trait_sd = sd(value, na.rm = T), n = n())

t1 <- gc.prm.stat %>%
    rename(rhizobia = strain) %>%
    filter(rhizobia %in% rhizobia_strains) %>%
    select(rhizobia, r, lag, maxOD) %>%
    pivot_longer(cols = -rhizobia, names_to = "trait", values_to = "trait_mean")
t2 <- gc.prm.stat %>%
    rename(rhizobia = strain) %>%
    filter(rhizobia %in% rhizobia_strains) %>%
    select(rhizobia, r=r.sem, lag=lag.sem, maxOD=maxOD.sem) %>%
    pivot_longer(cols = -rhizobia, names_to = "trait", values_to = "trait_sd")

gc_long_stat <- left_join(t1, t2) %>% mutate(n = 4) %>% mutate(rhizobia_site = str_sub(rhizobia, 1, 1))

trait_long_stat <- bind_rows(mutate(gc_long_stat, trait_type = "growth"), mutate(treatments_long_stat, trait_type = "mutualism")) %>%
    select(rhizobia_site, rhizobia, everything()) %>%
    arrange(rhizobia_site, rhizobia, trait)



# 1. plot growth trait mean by site ----
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
ggsave(paste0(folder_data, "temp/22-01-gc_trait_site.png"), p, width = 6, height = 5)

##
gc_h <- gc.prm.stat %>% filter(site == "H")
gc_l <- gc.prm.stat %>% filter(site == "L")

# lag
wilcox.test(gc_h$lag, gc_l$lag) # p=0.599
# r
wilcox.test(gc_h$r, gc_l$r) # p=0.4136
# maxOD
wilcox.test(gc_h$maxOD, gc_l$maxOD) # p=0.8518

# 2. pca of all growth traits ----


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



