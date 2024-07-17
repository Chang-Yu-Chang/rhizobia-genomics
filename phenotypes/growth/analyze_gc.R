#' This script analyzes the growth curve data

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gcs <- read_csv(paste0(folder_phenotypes, 'growth/gcs.csv'))
gc_summs <- read_csv(paste0(folder_phenotypes, 'growth/gc_summs.csv'))
gc_prms <- read_csv(paste0(folder_phenotypes, 'growth/gc_prms.csv'))
gc_prm_summs <- read_csv(paste0(folder_phenotypes, 'growth/gc_prm_summs.csv'))


# 0. clean up data ----
factorize_vars <- function (tb) tb %>% mutate(exp_id = factor(exp_id, isolates$exp_id), temperature = factor(temperature, c("25c", "30c", "35c", "40c")))
isolates <- isolates %>% mutate(exp_id = factor(exp_id, isolates$exp_id))
gc_summs <- factorize_vars(gc_summs)
gc_prms <- factorize_vars(gc_prms)
gc_summs <- factorize_vars(gc_summs)

isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD))


# 1. Plot all growth traits
p <- gcs %>%
    mutate(exp_id = factor(exp_id, isolates$exp_id)) %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = exp_id, group = well)) +
    facet_wrap(~temperature, nrow = 2) +
    scale_x_continuous(breaks = seq(0, 48, 12)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(color = "grey90")
    ) +
    guides() +
    labs(x = "time (hr)")

ggsave(paste0(folder_phenotypes, "growth/01-gc_raw.png"), p, width = 10, height = 8)

# 2. Growth curves by well
p <- gcs %>%
    mutate(exp_id = factor(exp_id, isolates$exp_id)) %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = exp_id, alpha = temperature), linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 48, 12)) +
    scale_alpha_manual(values = c("25c" = 0.1, "30c" = 0.4, "35c" = 0.7, "40c" = 1)) +
    #scale_color_manual(values = rep(brewer.pal(n = 9, name = "Set1"), 5)) +
    facet_grid(row ~ column) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(color = "grey90"),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs()

ggsave(paste0(folder_phenotypes, "growth/02-gc_raw_grid.png"), p, width = 20, height = 15)

# 3. Plot the pairwise comparison of the r
plot_dots <- function(tb_wide, trait1, trait2, axis_lower, axis_upper) {
    tb_wide %>%
        ggplot() +
        geom_smooth(aes(x = {{trait1}}, y = {{trait2}}), method = "lm")+
        geom_point(aes(x = {{trait1}}, y = {{trait2}}), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
        #scale_color_manual(values = brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)]) +
        #facet_wrap(~site_group, nrow = 2) +
        scale_x_continuous(limits = c(axis_lower, axis_upper)) +
        scale_y_continuous(limits = c(axis_lower, axis_upper)) +
        coord_equal() +
        theme_classic() +
        theme(
            panel.grid.major = element_line(color = "grey95")
        ) +
        guides()
    #    labs(x = "# of nodules", y = "shoot biomass (mg)")
}
p1 <- isolates_gc %>% plot_dots(r_30c, r_25c, -0.1, 1.4)
p2 <- isolates_gc %>% plot_dots(r_35c, r_25c, -0.1, 1.4)
p3 <- isolates_gc %>% plot_dots(r_40c, r_25c, -0.1, 1.4)
p4 <- isolates_gc %>% plot_dots(r_35c, r_30c, -0.1, 1.4)
p5 <- isolates_gc %>% plot_dots(r_40c, r_30c, -0.1, 1.4)
p6 <- isolates_gc %>% plot_dots(r_40c, r_35c, -0.1, 1.4)

p <- plot_grid(p1, p2, p3, NULL, p4, p5, NULL, NULL, p6,
    nrow = 3, align = "hv", axis = "tblr", scale = 0.95, labels = c("A", "B", "C", "", "D", "E", "", "", "F")
) + theme(plot.background = element_rect(color = NA, fill = "white"))
p
ggsave(paste0(folder_phenotypes, "growth/03-r_pairs.png"), p, width = 9, height = 9)

# 4. plot the lag pairs
range(c(isolates_gc$lag_25c, isolates_gc$lag_30c, isolates_gc$lag_35c ,isolates_gc$lag_40c), na.rm = T)
p1 <- isolates_gc %>% plot_dots(lag_30c, lag_25c, 5, 40)
p2 <- isolates_gc %>% plot_dots(lag_35c, lag_25c, 5, 40)
p3 <- isolates_gc %>% plot_dots(lag_40c, lag_25c, 5, 40)
p4 <- isolates_gc %>% plot_dots(lag_35c, lag_30c, 5, 40)
p5 <- isolates_gc %>% plot_dots(lag_40c, lag_30c, 5, 40)
p6 <- isolates_gc %>% plot_dots(lag_40c, lag_35c, 5, 40)

p <- plot_grid(p1, p2, p3, NULL, p4, p5, NULL, NULL, p6,
    nrow = 3, align = "hv", axis = "tblr", scale = 0.95, labels = c("A", "B", "C", "", "D", "E", "", "", "F")
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(paste0(folder_phenotypes, "growth/04-lag_pairs.png"), p, width = 9, height = 9)

# 4. plot the lag pairs
range(c(isolates_gc$maxOD_25c, isolates_gc$maxOD_30c, isolates_gc$maxOD_35c ,isolates_gc$maxOD_40c), na.rm = T)
p1 <- isolates_gc %>% plot_dots(maxOD_30c, maxOD_25c, 0.05, 0.45)
p2 <- isolates_gc %>% plot_dots(maxOD_35c, maxOD_25c, 0.05, 0.45)
p3 <- isolates_gc %>% plot_dots(maxOD_40c, maxOD_25c, 0.05, 0.45)
p4 <- isolates_gc %>% plot_dots(maxOD_35c, maxOD_30c, 0.05, 0.45)
p5 <- isolates_gc %>% plot_dots(maxOD_40c, maxOD_30c, 0.05, 0.45)
p6 <- isolates_gc %>% plot_dots(maxOD_40c, maxOD_35c, 0.05, 0.45)

p <- plot_grid(p1, p2, p3, NULL, p4, p5, NULL, NULL, p6,
    nrow = 3, align = "hv", axis = "tblr", scale = 0.95, labels = c("A", "B", "C", "", "D", "E", "", "", "F")
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(paste0(folder_phenotypes, "growth/05-maxOD_pairs.png"), p, width = 9, height = 9)

# Plot the two strains shared in mine and Linda's experiment
p <- gcs %>%
    filter(genome_id %in% c("g4", "g13")) %>%
    #filter(genome_id %in% paste0("g", c(4,5,6,9,11,13,16))) %>%
    mutate(group = paste0(temperature, well, genome_id)) %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = genome_id, group = group)) +
    #scale_color_manual(values = site_group_colors) +
    scale_color_manual(values = c(g4="#0C6291", g13="#BF4342")) +
    facet_grid(.~temperature) +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "phenotypes_analysis/growth/06-gc_two_strains.png"), p, width = 10, height = 3)

# Plot all strains used in Linda's experiment
p <- gcs %>%
    #filter(genome_id %in% c("g4", "g13")) %>%
    filter(genome_id %in% paste0("g", c(4,5,6,9,11,13,16))) %>%
    mutate(group = paste0(temperature, well, genome_id)) %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = site_group, group = group)) +
    scale_color_manual(values = site_group_colors) +
    #scale_color_manual(values = c(g4="#0C6291", g13="#BF4342")) +
    facet_grid(.~temperature) +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "phenotypes_analysis/growth/07-gc_lindas_strain.png"), p, width = 10, height = 3)

# All symbiontic strain's gc
p <- gcs %>%
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    filter(population == "VA") %>%
    #filter(genome_id %in% c("g4", "g13")) %>%
    mutate(group = paste0(temperature, well)) %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = site_group, group = group)) +
    scale_color_manual(values = site_group_colors) +
    #scale_color_manual(values = c(g4="#0C6291", g13="#BF4342")) +
    facet_grid(.~temperature) +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "phenotypes_analysis/growth/08-gc_symbiotic.png"), p, width = 10, height = 3)



