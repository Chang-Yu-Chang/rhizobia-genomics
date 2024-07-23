#' This script plot the growth traits

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))


gts <- read_csv(paste0(folder_data, "phenotypes/growth/gts.csv")) %>%
    clean_names() %>%
    filter(temperature != "40c") %>%
    left_join(isolates) %>%
    select(exp_id, r, lag, max_od, temperature) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, max_od)) %>%
    drop_na

# Plot the pairs
tt1 <- gts %>%
    left_join(isolates) %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    filter(population == "VA")
tt1_summ <- tt1 %>% group_by(site_group) %>% summarize(mean = mean(r_30c), sem = sd(r_30c)/sqrt(n()))
p1 <- tt1 %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = r_30c, fill = site_group)) +
    #geom_segment(data = tt1_summ, aes(x = site_group, xend = site_group, y = mean-sem*2, yend = mean+sem*2), linewidth = 1, color = "black") +
    geom_jitter(aes(x = site_group, y = r_30c, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    facet_wrap(population ~., scale = "free_x", nrow = 2) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA)
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "growth rate at 30C (1/hr)")

tt2 <- gts %>%
    left_join(isolates) %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    filter(population == "PA")
tt2_summ <- tt2 %>% group_by(site_group) %>% summarize(mean = mean(r_30c), sem = sd(r_30c)/sqrt(n()))

p2 <- gts %>%
    left_join(isolates) %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    filter(population == "PA") %>%
    ggplot() +
    geom_boxplot(aes(x = site_group, y = r_30c, fill = site_group)) +
    #geom_segment(data = tt2_summ, aes(x = site_group, xend = site_group, y = mean-sem*2, yend = mean+sem*2), linewidth = 1, color = "black") +
    geom_jitter(aes(x = site_group, y = r_30c, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
    scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
    scale_color_manual(values = alpha(site_group_colors, 1)) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    facet_wrap(population ~., scale = "free_x", nrow = 2) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(color = "grey10", fill = NA)
    ) +
    guides(fill = "none", color = "none") +
    labs(y = "growth rate at 30C (1/hr)")

# Plot the group
gts1 <- gts[gts$exp_id %in% isolates$exp_id[isolates$population == "VA"],]
gts2 <- gts[gts$exp_id %in% isolates$exp_id[isolates$population == "PA"],]
pca_results1 <- prcomp(gts1[,-1], scale. = TRUE)
pca_results2 <- prcomp(gts2[,-1], scale. = TRUE)
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

pcs <- bind_rows(
    as_tibble(pca_results1$x) %>% mutate(exp_id = gts1$exp_id),
    as_tibble(pca_results2$x) %>% mutate(exp_id = gts2$exp_id)
) %>%
    left_join(isolates)

p3 <- pcs %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    filter(population == "VA") %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,4,2)) +
    scale_y_continuous(breaks = seq(-2,2,2)) +
    facet_wrap(population ~., nrow = 2, scales = "free_x") +
    theme_minimal() +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1)
    ) +
    guides(color = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results1)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results1)[2], "%)"))

p4 <-  pcs %>%
    mutate(population = factor(population, c("VA", "PA"))) %>%
    filter(population == "PA") %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2) +
    geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(breaks = seq(-4,4,2)) +
    scale_y_continuous(breaks = seq(-2,2,2)) +
    facet_wrap(population ~., nrow = 2, scales = "free_x") +
    theme_minimal() +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey10", fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1)
    ) +
    guides(color = "none") +
    labs(x = paste0("PC1 (", get_pcvar(pca_results2)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results2)[2], "%)"))


p <- plot_grid(p1, p3, p2, p4, align = T, axis = "tb", scale = 0.95, rel_widths = c(1, 2))

ggsave(here::here("forposter/growth.pdf"), p, width = 4, height = 6)
