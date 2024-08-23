#' This script

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for permanova
source(here::here("metadata.R"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>%
    slice(1:32)

set.seed(1)
detect_sig <- function (pv) {
    if (pv > 0.05) {
        return("n.s.")
    } else if (pv > 0.01) {
        return("*")
    } else if (pv > 0.001) {
        return("**")
    } else if (pv < 0.001) {
        return("***")
    }
}

# Growth rate ----
gts <- read_csv(paste0(folder_data, "phenotypes/growth/gts.csv")) %>%
    clean_names() %>%
    filter(temperature != "40c") %>%
    #left_join(isolates) %>%
    select(exp_id, r, lag, max_od, temperature) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, max_od)) %>%
    drop_na %>%
    left_join(isolates)

# r_30c ~ elevation w/ N
tb <- gts %>%
    left_join(isolates) %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    filter(population == "VA")
mod <- lmer(r_30c ~ site_group + (1|site), data = tb)
Anova(mod, type = 3)
vmax <- max(tb$r_30c)*.95

p1 <- tb %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = r_30c, fill = site_group)) +
    geom_jitter(aes(x = site_group, y = r_30c, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
    annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
    annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    #scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_text(size = 10),
        panel.border = element_rect(color = "grey10", fill = NA),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "growth rate at 30C (1/hr)")


# r_30c ~ elevation w/ N
tb <- gts %>%
    left_join(isolates) %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    filter(population == "PA")
mod <- lmer(r_30c ~ site_group + (1|site), data = tb)
Anova(mod, type = 3)
vmax <- max(tb$r_30c)*.95

p2 <- tb %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = r_30c, fill = site_group)) +
    geom_jitter(aes(x = site_group, y = r_30c, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
    annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
    annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_text(size = 10),
        panel.border = element_rect(color = "grey10", fill = NA),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "growth rate at 30C (1/hr)")

# PCA
gts1 <- gts %>%
    filter(population == "VA") %>%
    select(site_group, starts_with(c("r", "lag", "max"))) %>%
    drop_na()
gts2 <- gts %>%
    filter(population == "PA") %>%
    select(site_group, starts_with(c("r", "lag", "max"))) %>%
    drop_na()
pca_results1 <- prcomp(gts1[,-1], scale. = TRUE)
pca_results2 <- prcomp(gts2[,-1], scale. = TRUE)
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

pcs <- bind_rows(
    as_tibble(pca_results1$x) %>% mutate(site_group = gts1$site_group),
    as_tibble(pca_results2$x) %>% mutate(site_group = gts2$site_group)
) %>%
    left_join(distinct(isolates, population, site_group)) %>%
    mutate(population = factor(population, c("VA", "PA")))

# traits ~ elevation
pcsi <- pcs %>% filter(population == "VA")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- adonis2(dm ~ site_group, data = pcsi)
mod

p3 <- pcsi %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-8,8,2)) +
    scale_y_continuous(breaks = seq(-8,8,2)) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        plot.background = element_blank()
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results1)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results1)[2], "%)"))

# traits ~ elevation
pcsi <- pcs %>% filter(population == "PA")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- adonis2(dm ~ site_group, data = pcsi)
mod

p4 <- pcsi %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-8,8,2)) +
    scale_y_continuous(breaks = seq(-6,6,2)) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        plot.background = element_blank()
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results2)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results2)[2], "%)"))

p_growth <- plot_grid(p1, p3, p2, p4, nrow = 1, align = T, axis = "tb", scale = 0.85, rel_widths = c(1, 2), labels = c("A", "", "B", ""))

# Lupulina experiment ----
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))

# shoot_biomass ~ elevation w/o N
tb <- plants %>%
    filter(exp_plant == "lupulina", exp_id != "control", population == "VA") %>%
    drop_na(shoot_biomass_mg)
mod <- lmer(shoot_biomass_mg ~ site_group + (1|genome_id), data = tb)
Anova(mod, type = 3)
vmax <- max(tb$shoot_biomass_mg)*.95

p1 <- tb %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_biomass_mg, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_biomass_mg, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
    annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot biomass (mg)")

# shoot_biomass ~ elevation w/o N
tb <- plants %>%
    filter(exp_plant == "lupulina", exp_id != "control", population == "PA") %>%
    drop_na(shoot_biomass_mg)
mod <- lmer(shoot_biomass_mg ~ site_group + (1|genome_id), data = tb)
Anova(mod, type = 3)
vmax <- max(tb$shoot_biomass_mg)*.95

p2 <- tb %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_biomass_mg, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_biomass_mg, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
    annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot biomass (mg)")

# plot the PCA for lupulina
plants1 <- plants %>%
    filter(exp_plant == "lupulina", exp_id != "control", population == "VA") %>%
    select(site_group, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
    drop_na()
plants2 <- plants %>%
    filter(exp_plant == "lupulina", exp_id != "control", population == "PA") %>%
    select(site_group, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
    drop_na()

pca_results1 <- prcomp(plants1[,-1], scale. = TRUE)
pca_results2 <- prcomp(plants2[,-1], scale. = TRUE)
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

pcs <- bind_rows(
    as_tibble(pca_results1$x) %>% mutate(site_group = plants1$site_group),
    as_tibble(pca_results2$x) %>% mutate(site_group = plants2$site_group)
) %>%
    left_join(distinct(isolates, population, site_group))

# traits ~ site_group
pcsi <- pcs %>% filter(population == "VA")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- adonis2(dm ~ site_group, data = pcsi)
mod

p3 <- pcsi %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2, alpha = .5) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(plants1)), hjust = 1.1, vjust = 1.1) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,8,2)) +
    scale_y_continuous(breaks = seq(-2,2,1)) +
    theme_bw() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        plot.background = element_blank()
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results1)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results1)[2], "%)"))

# traits ~ site_group
pcsi <- pcs %>% filter(population == "PA")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- adonis2(dm ~ site_group, data = pcsi)
mod

p4 <- pcsi %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2, alpha = .5) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(plants2)), hjust = 1.1, vjust = 1.1) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,6,2)) +
    scale_y_continuous(breaks = seq(-2,2,1)) +
    theme_bw() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        plot.background = element_blank()
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results2)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results2)[2], "%)"))

p_lupulina <- plot_grid(p1, p3, p2, p4, nrow = 1, align = T, axis = "tb", scale = 0.85, rel_widths = c(1, 2), labels = c("C", "", "D", ""))


# Sativa experiment ----
# shoot_height ~ elevation w/o N
tb <- plants %>%
    filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", population == "VA") %>%
    drop_na(shoot_height)
mod <- lmer(shoot_height ~ site_group + (1|genome_id), data = tb)
Anova(mod, type = 3)
vmax <- max(tb$shoot_height)*.95/10

p1 <- tb %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_height/10, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_height/10, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
    annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot height (cm)")


# shoot_height ~ elevation w/o N
tb <- plants %>%
    filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", population == "PA") %>%
    drop_na(shoot_height)
mod <- lmer(shoot_height ~ site_group + (1|genome_id), data = tb)
Anova(mod, type = 3)
vmax <- max(tb$shoot_height)*.95

p2 <- tb %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = shoot_height, fill = site_group), outliers = F) +
    geom_jitter(aes(x = site_group, y = shoot_height, color = site_group), size = 2, width = .1, shape = 21, stroke = 1, alpha = .5) +
    annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
    annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA),
        plot.background = element_blank()
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "shoot height (cm)")


# plot the PCA for sativa
plants1 <- plants %>%
    filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", population == "VA") %>%
    select(site_group, shoot_height, nodule_number, longest_petiole_length, leaf_number, leaf_color) %>%
    drop_na()
plants2 <- plants %>%
    filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", population == "PA") %>%
    select(site_group, shoot_height, nodule_number, leaf_number, leaf_color, lateral_root_number, longest_lateral_root_length) %>%
    drop_na()

pca_results1 <- prcomp(plants1[,-1], scale. = TRUE)
pca_results2 <- prcomp(plants2[,-1], scale. = TRUE)
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

# traits ~ site_group
pcsi <- as_tibble(pca_results1$x) %>%
    mutate(site_group = plants1$site_group) %>%
    left_join(distinct(isolates, population, site_group))
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- adonis2(dm ~ site_group, data = pcsi)
mod

p3 <- pcsi %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2, alpha = .5) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(plants1)), hjust = 1.1, vjust = 1.1) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,8,2)) +
    scale_y_continuous(breaks = seq(-4,10,2)) +
    theme_bw() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        plot.background = element_blank()
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results1)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results1)[2], "%)"))

# traits ~ site_group
pcsi <- as_tibble(pca_results2$x) %>%
    mutate(site_group = plants2$site_group) %>%
    left_join(distinct(isolates, population, site_group))
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- adonis2(dm ~ site_group, data = pcsi)
mod

p4 <- pcsi %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2, alpha = .5) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(plants2)), hjust = 1.1, vjust = 1.1) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_fill_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,6,2)) +
    scale_y_continuous(breaks = seq(-4,4,2)) +
    theme_bw() +
    theme(
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        plot.background = element_blank()
    ) +
    guides(color = "none", fill = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results2)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results2)[2], "%)"))


p_sativa <- plot_grid(p1, p3, p2, p4, nrow = 1, align = T, axis = "tb", scale = 0.85, rel_widths = c(1, 2), labels = c("E", "", "F", ""))



# Combine figures ----

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig2.png"), scale = 1) +
    draw_plot(p_growth, width = .6, height = .32, x = .35, y = .59) +
    draw_plot(p_lupulina, width = .6, height = .32, x = .35, y = .29) +
    draw_plot(p_sativa, width = .6, height = .32, x = .35, y = -0.01) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig2.png"), p, width = 14, height = 8)

