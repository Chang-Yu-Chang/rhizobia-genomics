#' This script plot the symbiosis traits

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))

set.seed(1)
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))


# Plot the pairs
p1 <- plants %>%
    # group_by(exp_plant, population, site_group, exp_id, exp_nitrogen) %>%
    # summarize(shoot_biomass_mg = mean(shoot_biomass_mg, na.rm = T))
    filter(exp_plant == "sativa") %>%
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    drop_na(shoot_height) %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_height/10, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_height/10, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA)
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot height (cm)")


# Plot the pairs
p2 <- plants %>%
    filter(exp_plant == "sativa") %>%
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    drop_na(shoot_height) %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_height, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_height, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA)
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot height (cm)")


# plot the PCA for sativa
plants1 <- plants %>%
    filter(exp_plant == "sativa") %>%
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    #mutate(id = 1:n()) %>%
    select(site_group, shoot_height, nodule_number, longest_petiole_length, leaf_number, leaf_color) %>%
    drop_na()
plants2 <- plants %>%
    filter(exp_plant == "sativa") %>%
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    select(site_group, shoot_height, nodule_number, leaf_number, leaf_color, lateral_root_number, longest_lateral_root_length) %>%
    drop_na()
dim(plants1)
dim(plants2)

pca_results1 <- prcomp(plants1[,-1], scale. = TRUE)
pca_results2 <- prcomp(plants2[,-1], scale. = TRUE)
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

pcs <- bind_rows(
    as_tibble(pca_results1$x) %>% mutate(site_group = plants1$site_group),
    as_tibble(pca_results2$x) %>% mutate(site_group = plants2$site_group)
)

p3 <- pcs %>%
    filter(site_group %in% c("high elevation", "low elevation")) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2, alpha = .5) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(plants1)), hjust = 1.1, vjust = 1.1) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,8,2)) +
    scale_y_continuous(breaks = seq(-4,10,2)) +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1)
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results1)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results1)[2], "%)"))

p4 <- pcs %>%
    filter(site_group %in% c("urban", "suburban")) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2, alpha = .5) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(plants2)), hjust = 1.1, vjust = 1.1) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,6,2)) +
    scale_y_continuous(breaks = seq(-4,4,2)) +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1)
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results2)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results2)[2], "%)"))


p <- plot_grid(p1, p3, p2, p4, align = T, axis = "tb", scale = 0.95, rel_widths = c(1, 2))

ggsave(here::here("forposter/plants_sat.pdf"), p, width = 4, height = 6)

#
sativa <- plants %>%
    filter(exp_plant == "sativa") %>%
    filter(exp_id != "control") %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    drop_na(shoot_height)

sam <- sativa %>%
    group_by(population, site_group) %>%
    count()

p <- sativa %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_height, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_height, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    geom_text(data = sam, aes(x = site_group, label = paste0("N=", n)), y = Inf, vjust = 1.2, size = 3) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    scale_y_continuous(expand = c(.1,.1)) +
    facet_wrap(.~population, scales = "free") +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(color = "grey10", fill = NA),
        strip.text = element_blank(),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot height (cm)")

 ggsave(here::here("forposter/plants_sat.pdf"), p, width = 4, height = 4)
