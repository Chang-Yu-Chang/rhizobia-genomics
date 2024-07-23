#' This script plot the symbiosis traits

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))


plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))

# Plot the pairs
p1 <- plants %>%
    # group_by(exp_plant, population, site_group, exp_id, exp_nitrogen) %>%
    # summarize(shoot_biomass_mg = mean(shoot_biomass_mg, na.rm = T))
    filter(exp_plant == "lupulina") %>%
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    drop_na(shoot_biomass_mg) %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_biomass_mg, fill = site_group)) +
    geom_jitter(aes(x = site_group, y = shoot_biomass_mg, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    #scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    facet_wrap(population ~., scale = "free_x", nrow = 2) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA)
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot biomass (mg)")

p2 <- plants %>%
    # group_by(exp_plant, population, site_group, exp_id, exp_nitrogen) %>%
    # summarize(shoot_biomass_mg = mean(shoot_biomass_mg, na.rm = T))
    filter(exp_plant == "lupulina") %>%
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    drop_na(shoot_biomass_mg) %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_biomass_mg, fill = site_group)) +
    geom_jitter(aes(x = site_group, y = shoot_biomass_mg, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    #scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    facet_wrap(population ~., scale = "free_x", nrow = 2) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA)
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot biomass (mg)")

# plot the PCA for lupulina
plants1 <- plants %>%
    filter(exp_plant == "lupulina") %>%
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    #mutate(id = 1:n()) %>%
    select(site_group, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
    drop_na()
plants2 <- plants %>%
    filter(exp_plant == "lupulina") %>%
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    select(site_group, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
    drop_na()


# gts1 <- plants[plants$exp_id %in% isolates$exp_id[isolates$population == "VA"],]
# gts2 <- gts[gts$exp_id %in% isolates$exp_id[isolates$population == "PA"],]
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
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,8,2)) +
    scale_y_continuous(breaks = seq(-2,2,2)) +
    #facet_wrap(population ~., nrow = 2, scales = "free_x") +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1)
    ) +
    guides(color = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results1)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results1)[2], "%)"))

p4 <- pcs %>%
    filter(site_group %in% c("urban", "suburban")) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,6,2)) +
    scale_y_continuous(breaks = seq(-2,2,2)) +
    #facet_wrap(population ~., nrow = 2, scales = "free_x") +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1)
    ) +
    guides(color = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results1)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results1)[2], "%)"))


p <- plot_grid(p1, p3, p2, p4, align = T, axis = "tb", scale = 0.95, rel_widths = c(1, 2))

ggsave(here::here("forposter/plants.pdf"), p, width = 4, height = 6)

