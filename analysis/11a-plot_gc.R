#' This script fits the growth curve data

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"), show_col_types = F)
gc <- read_csv(paste0(folder_data, 'temp/11-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/11-gc_summ.csv'), show_col_types = F)
gc_prm <- read_csv(paste0(folder_data, 'temp/11-gc_prm.csv'), show_col_types = F)
gc_prm_summ <- read_csv(paste0(folder_data, 'temp/11-gc_prm_summ.csv'), show_col_types = F)

# 0. clean up data ----
isolates_mapping <- isolates_mapping %>% mutate(exp_id = factor(exp_id, isolates_mapping$exp_id))
gc_summ <- gc_summ %>% mutate(exp_id = factor(exp_id, isolates_mapping$exp_id))
gc_prm <- gc_prm %>% mutate(exp_id = factor(exp_id, isolates_mapping$exp_id))
gc_prm_summ <- gc_prm_summ %>% mutate(exp_id = factor(exp_id, isolates_mapping$exp_id))

# 1. Raw data ----
p <- gc %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = exp_id, group = well)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/11a-01-gc_raw.png"), p, width = 8, height = 5)


# 2. Growth curves by well ----
p <- gc %>%
    ggplot() +
    geom_line(aes(x = t, y = abs, color = exp_id), linewidth = 1) +
    theme_light() +
    facet_grid(row ~ column) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/11a-02-gc_raw_grid.png"), p, width = 20, height = 15)

# 3. Average OD by exp_id ----
p <- gc_summ %>%
    left_join(isolates_mapping) %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs, color = rhizobia_site)) +
    geom_ribbon(aes(x = t, ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs, fill = rhizobia_site), alpha = 0.2) +
    facet_wrap(~exp_id, ncol = 4) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA)
    ) +
    #guides(color = "none", fill = "none") +
    labs(x = "time (hrs)", y = expression(OD[600]))
ggsave(paste0(folder_data, "temp/11a-03-gc_mean.png"), p, width = 8, height = 12)

# 4. Average OD by exp_id, overlayed ----
p <- gc_summ %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(group = exp_id)) +
    geom_line(aes(x = t, y = mean_abs, color = rhizobia_site)) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "right"
    ) +
    labs(x = "time (hrs)", y = expression(OD[600]))
ggsave(paste0(folder_data, "temp/11a-04-gc_mean_overlay.png"), p, width = 6, height = 4)


# 5. plot all trait mean by strain_site_group ----
p1 <- gc_prm_summ %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(x = rhizobia_site, y = lag, fill = rhizobia_site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0, height = 0, size = 2, stroke = 1) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = "lag time (hr)")

p2 <- gc_prm_summ %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(x = rhizobia_site, y = r, fill = rhizobia_site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0, height = 0, size = 2, stroke = 1) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = expression(growth~rate(h^-1)))

p3 <- gc_prm_summ %>%
    left_join(isolates_mapping) %>%
    ggplot(aes(x = rhizobia_site, y = maxOD, fill = rhizobia_site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0, height = 0, size = 2, stroke = 1) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        legend.position = "right"
    ) +
    guides() +
    labs(x = "", y = expression(paste("max", "[", OD[600], "]")))

p <- plot_grid(p1, p2, p3, nrow = 1, axis = "tbrl", align = "h", rel_widths = c(1,1,1.6))
ggsave(paste0(folder_data, "temp/11a-05-gc_trait_site.png"), p, width = 8, height = 4)


# 6. check the taxonomy ----
gc_prm_summ %>%
    left_join(isolates_mapping) %>%
    filter(rhizobia_site %in% c("high-elevation", "low-elevation")) %>%
    arrange(desc(maxOD)) %>%
    select(exp_id, maxOD, rhizobia_site, genome_id)






#
# ##
# gc_h <- gc_prm_summ %>% filter(strain_site_group == "H")
# gc_l <- gc_prm_summ %>% filter(strain_site_group == "L")
#
# # lag
# wilcox.test(gc_h$lag, gc_l$lag) # p=0.599
# # r
# wilcox.test(gc_h$r, gc_l$r) # p=0.4136
# # maxOD
# wilcox.test(gc_h$maxOD, gc_l$maxOD) # p=0.8518
#
#
# # 8. plot all trait raw by strain_site_group ----
# p1 <- gc_prm %>%
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
# p2 <- gc_prm %>%
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
# p3 <- gc_prm %>%
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
# ggsave(paste0(folder_data, "temp/11a-08-gc_trait_site.png"), p, width = 6, height = 5)


