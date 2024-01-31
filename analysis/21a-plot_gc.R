#' This script plots the growth curve data

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"))
gcs <- read_csv(paste0(folder_data, 'temp/21-gcs.csv'))
gc_summs <- read_csv(paste0(folder_data, 'temp/21-gc_summs.csv'))
gc_prms <- read_csv(paste0(folder_data, 'temp/21-gc_prms.csv'))
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'))


# 0. clean up data
factorize_vars <- function (tb) tb %>% mutate(exp_id = factor(exp_id, isolates_mapping$exp_id), temperature = factor(temperature, c("25c", "30c", "35c", "40c")))
isolates_mapping <- isolates_mapping %>% mutate(exp_id = factor(exp_id, isolates_mapping$exp_id))
gc_summs <- factorize_vars(gc_summs)
gc_prms <- factorize_vars(gc_prms)
gc_summs <- factorize_vars(gc_summs)

isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD))



# 1. Plot all growth traits
p <- gcs %>%
    mutate(exp_id = factor(exp_id, isolates_mapping$exp_id)) %>%
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

ggsave(paste0(folder_data, "temp/21a-01-gc_raw.png"), p, width = 10, height = 8)

# 2. Growth curves by well 
p <- gcs %>%
    mutate(exp_id = factor(exp_id, isolates_mapping$exp_id)) %>%
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

ggsave(paste0(folder_data, "temp/21a-02-gc_raw_grid.png"), p, width = 20, height = 15)

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
ggsave(paste0(folder_data, "temp/21a-03-r_pairs.png"), p, width = 9, height = 9)

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

ggsave(paste0(folder_data, "temp/21a-04-lag_pairs.png"), p, width = 9, height = 9)

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

ggsave(paste0(folder_data, "temp/21a-05-maxOD_pairs.png"), p, width = 9, height = 9)



if (FALSE) {

# 3. Average OD by exp_id 
p <- gc_summs %>%
    left_join(isolates_mapping) %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs, color = rhizobia_site, alpha = temperature)) +
    scale_x_continuous(breaks = seq(0, 48, 12)) +
    scale_alpha_manual(values = c("25c" = 0.1, "30c" = 0.4, "35c" = 0.7, "40c" = 1)) +
    #geom_ribbon(aes(x = t, ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs, fill = rhizobia_site), alpha = 0.2) +
    facet_wrap(~exp_id, ncol = 4) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        panel.grid.major.x = element_line(color = "grey90")
    ) +
    #guides(color = "none", fill = "none") +
    labs(x = "time (hrs)", y = expression(OD[600]))
ggsave(paste0(folder_data, "temp/21a-03-gc_mean.png"), p, width = 8, height = 12)

# 4. Average OD by exp_id, overlayed 
p <- gc_summs %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(group = exp_id)) +
    geom_line(aes(x = t, y = mean_abs, color = rhizobia_site, group = exp_id)) +
    scale_color_manual(values = rhizobia_site_colors) +
    facet_grid(.~temperature) +
    #scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        panel.grid.major.x = element_line(color = "grey90"),
        legend.position = "top"
    ) +
    labs(x = "time (hrs)", y = expression(OD[600]))
ggsave(paste0(folder_data, "temp/21a-04-gc_mean_overlay.png"), p, width = 12, height = 4)


# 5. plot all trait mean by strain_site_group 
p1 <- gc_prm_summs %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(x = rhizobia_site, y = lag, fill = rhizobia_site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0, height = 0, size = 2, stroke = 1) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    facet_grid(.~temperature) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = "lag time (hr)")

p2 <- gc_prm_summs %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(x = rhizobia_site, y = r, fill = rhizobia_site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0, height = 0, size = 2, stroke = 1) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    facet_grid(.~temperature) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = expression(growth~rate(h^-1)))

p3 <- gc_prm_summs %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(x = rhizobia_site, y = maxOD, fill = rhizobia_site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0, height = 0, size = 2, stroke = 1) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    facet_grid(.~temperature) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        legend.position = "right"
    ) +
    guides() +
    labs(x = "", y = expression(paste("max", "[", OD[600], "]")))

p <- plot_grid(p1, p2, p3, nrow = 3, axis = "tbrl", align = "v")
ggsave(paste0(folder_data, "temp/21a-05-gc_trait_site.png"), p, width = 8, height = 12)

# 6. reaction norm 
p1 <- gc_prm_summs %>%
    left_join(isolates_mapping) %>%
    ggplot() +
    geom_line(aes(x = temperature, y = lag, group = exp_id, color = rhizobia_site)) +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        legend.position = "none"
    ) +
    guides() +
    labs(x = "", y = "lag time (hr)")

p2 <- gc_prm_summs %>%
    left_join(isolates_mapping) %>%
    ggplot() +
    geom_line(aes(x = temperature, y = r, group = exp_id, color = rhizobia_site)) +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1)
    ) +
    guides() +
    labs(x = "", y = expression(growth~rate(h^-1)))

p3 <- gc_prm_summs %>%
    left_join(isolates_mapping) %>%
    ggplot() +
    geom_line(aes(x = temperature, y = maxOD, group = exp_id, color = rhizobia_site)) +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        legend.position = "none"
    ) +
    guides() +
    labs(x = "", y = expression(paste("max", "[", OD[600], "]")))

p <- plot_grid(p1, p2, p3, nrow = 3, axis = "tbrl", align = "v")
ggsave(paste0(folder_data, "temp/21a-06-reaction_norm.png"), p, width = 6, height = 12)


# 6. check the taxonomy 
gc_prm_summs %>%
    left_join(isolates_mapping) %>%
    filter(rhizobia_site %in% c("high-elevation", "low-elevation")) %>%
    arrange(desc(maxOD)) %>%
    select(exp_id, maxOD, rhizobia_site, genome_id)






#
# ##
# gc_h <- gc_prm_summs %>% filter(strain_site_group == "H")
# gc_l <- gc_prm_summs %>% filter(strain_site_group == "L")
#
# # lag
# wilcox.test(gc_h$lag, gc_l$lag) # p=0.599
# # r
# wilcox.test(gc_h$r, gc_l$r) # p=0.4136
# # maxOD
# wilcox.test(gc_h$maxOD, gc_l$maxOD) # p=0.8518
#
#
# # 8. plot all trait raw by strain_site_group 
# p1 <- gc_prms %>%
#     ggplot(aes(x = strain_site_group, y = lag, fill = strain_site_group)) +
#     geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
#     geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
#     scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
#     coord_flip() +
#     theme_classic() +
#     theme(
#         panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
#         axis.text.x = element_blank()
#     ) +
#     guides(fill = "none") +
#     labs(x = "", y = "lag time (hr)")
#
# p2 <- gc_prms %>%
#     ggplot(aes(x = strain_site_group, y = r, fill = strain_site_group)) +
#     geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
#     geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
#     scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
#     coord_flip() +
#     theme_classic() +
#     theme(
#         panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
#         axis.text.x = element_blank(),
#         legend.position = "right"
#     ) +
#     guides() +
#     labs(x = "", y = expression(growth~rate(h^-1)))
#
# p3 <- gc_prms %>%
#     ggplot(aes(x = strain_site_group, y = maxOD, fill = strain_site_group)) +
#     geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
#     geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
#     scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
#     coord_flip() +
#     theme_classic() +
#     theme(
#         panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
#         axis.text.x = element_blank()
#     ) +
#     guides(fill = "none") +
#     labs(x = "", y = expression(paste("max", "[", OD[600], "]")))
#
# p <- plot_grid(p1, p2, p3, ncol = 1, axis = "tbrl", align = "hv")
# ggsave(paste0(folder_data, "temp/21a-08-gc_trait_site.png"), p, width = 6, height = 5)


}

range(c(isolates_gc$r_25c,isolates_gc$r_30c,isolates_gc$r_35c,isolates_gc$r_40c))
