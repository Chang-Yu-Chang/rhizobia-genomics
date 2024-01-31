#' This script plots the trait correlation 

renv::load()
library(tidyverse)
library(janitor)
library(ggsci)
library(cowplot)
source(here::here("analysis/00-metadata.R"))


tb_traits <- read_csv(paste0(folder_data, "temp/24-tb_traits.csv"))
tb_traits_mean <- read_csv(paste0(folder_data, "temp/24-tb_traits_mean.csv"))
# tb_traits_dif <- read_csv(paste0(folder_data, "temp/24-tb_traits_dif.csv"))
# tb_reps_exp1 <- read_csv(paste0(folder_data, "temp/24-tb_reps_exp1.csv"))
# tb_reps_exp2 <- read_csv(paste0(folder_data, "temp/24-tb_reps_exp2.csv"))

# 
length(unique(tb_traits$trait1)) # 12 growth traits
length(unique(tb_traits$trait2)) # 37 symbiosis traits
# 12*37=444

# 1. plot the correaltion of one pair of trait  
name_trait1 <- tb_reps_exp1$trait1[1]
name_trait2 <- tb_reps_exp1$trait2[1]
colors_pair <- c(across = "grey", within = "maroon")
tb_reps_count <- tb_reps_exp1 %>%
    group_by(pair, trait1_value, trait2_value) %>%
    count()

# 1.1 trait correlation for only within pair strains
p1 <- tb_reps_count %>%
    filter(pair == "within") %>%
    ggplot() +
    geom_point(aes(x = trait1_value, y = trait2_value, color = pair, size = n), alpha = 0.3, stroke = 0) +
    scale_color_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = name_trait1, y = name_trait2)

# 1.2 trait correlation for both within-pair and across-pair
p2 <- tb_reps_count %>%
    ggplot() +
    geom_point(aes(x = trait1_value, y = trait2_value, color = pair, size = n), alpha = 0.3, stroke = 0) +
    scale_color_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = name_trait1, y = name_trait2)


# 1.3 line plot correlation coefficients of 100 bootstraps
tb_traits_long <- tb_traits %>%
    filter(trait1 == name_trait1[1], trait2 == name_trait2[1]) %>%
    pivot_longer(cols = c(cor_within, cor_across), names_to = "pair", names_prefix = "cor_")

p3 <- tb_traits_long %>%
    ggplot() +
    geom_point(aes(x = pair, y = value, group = bootstrap, color = pair), size = 3) +
    geom_line(aes(x = pair, y = value, group = bootstrap), linewidth = 0.5) +
    scale_color_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = "", y = "Spearman's rho")

# 1.4 histogram of plot the pairwise comparision
p4 <- tb_traits_long %>%
    ggplot() +
    geom_histogram(aes(x = value, fill = pair), position = "identity", alpha = 0.3) +
    geom_vline(xintercept = 0) +
    coord_flip() +
    scale_fill_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(fill = guide_legend(title = NULL)) +
    labs(x = "Pearson's r", y = "count")


# 
p <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4],
        scale = 0.95, align = "h", axis = "tb") +
        paint_white_background()


ggsave(paste0(folder_data, "temp/24a-01-correlation.png"), p, width = 10, height = 10)

tb_traits_dif %>%
    filter(trait1 == name_trait1[1], trait2 == name_trait2[1]) 

# 2. Plot another pair
name_trait1 <- tb_reps_exp2$trait1[1]
name_trait2 <- tb_reps_exp2$trait2[1]
colors_pair <- c(across = "grey", within = "maroon")
tb_reps_count <- tb_reps_exp2 %>%
    group_by(pair, trait1_value, trait2_value) %>%
    count()

# 2.1 trait correlation for only within pair strains
p1 <- tb_reps_count %>%
    filter(pair == "within") %>%
    ggplot() +
    geom_point(aes(x = trait1_value, y = trait2_value, color = pair, size = n), alpha = 0.3, stroke = 0) +
    scale_color_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = name_trait1, y = name_trait2)

# 2.2 trait correlation for both within-pair and across-pair
p2 <- tb_reps_count %>%
    ggplot() +
    geom_point(aes(x = trait1_value, y = trait2_value, color = pair, size = n), alpha = 0.3, stroke = 0) +
    scale_color_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = name_trait1, y = name_trait2)


# 2.3 line plot correlation coefficients of 100 bootstraps
tb_traits_long <- tb_traits %>%
    filter(trait1 == name_trait1[1], trait2 == name_trait2[1]) %>%
    pivot_longer(cols = c(cor_within, cor_across), names_to = "pair", names_prefix = "cor_")
p3 <- tb_traits_long %>%
    ggplot() +
    geom_point(aes(x = pair, y = value, group = bootstrap, color = pair), size = 3) +
    geom_line(aes(x = pair, y = value, group = bootstrap), linewidth = 0.5) +
    scale_color_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = "", y = "Spearman's rho")

# 2.4 histogram of plot the pairwise comparision
p4 <- tb_traits_long %>%
    ggplot() +
    geom_histogram(aes(x = value, fill = pair), position = "identity", alpha = 0.3) +
    geom_vline(xintercept = 0) +
    coord_flip() +
    scale_fill_manual(values = colors_pair) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(fill = guide_legend(title = NULL)) +
    labs(x = "Pearson's r", y = "count")


# 
p <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4],
        scale = 0.95, align = "h", axis = "tb") +
        paint_white_background()


ggsave(paste0(folder_data, "temp/24a-02-correlation.png"), p, width = 10, height = 10)

tb_traits_dif %>%
    filter(trait1 == name_trait1[1], trait2 == name_trait2[1]) 

# 3. Plot mean value
p1 <- tb_traits_mean %>%
    ggplot() +
    geom_point(aes(x = r_30c, y = dry_weight_mg), size = 3, shape = 21, stroke = 1) +
    geom_segment(aes(x = r_30c, xend = r_30c, y = dry_weight_mg - dry_weight_mg_se, yend = dry_weight_mg + dry_weight_mg_se)) +
    geom_segment(aes(x = r_30c - r_se_30c, xend = r_30c + r_se_30c, y = dry_weight_mg, yend = dry_weight_mg)) +
    theme_classic()

p2 <- tb_traits_mean %>%
    ggplot() +
    geom_point(aes(x = r_30c, y = nodule_number), size = 3, shape = 21, stroke = 1) +
    geom_segment(aes(x = r_30c, xend = r_30c, y = nodule_number - nodule_number_se, yend = nodule_number + nodule_number_se)) +
    geom_segment(aes(x = r_30c - r_se_30c, xend = r_30c + r_se_30c, y = nodule_number, yend = nodule_number)) +
    theme_classic()

p <- plot_grid(p1, p2, nrow = 1, scale = 0.95) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(paste0(folder_data, "temp/24a-03-trait_mean.png"), p, width = 8, height = 4)

# # 3. Plot the trait pairs that are correlated 
# tb_traits_dif %>%
#     filter(n >= 95)  %>%
    
    # ggplot() +
    # geom_histogram(aes(x = value, fill = name), position = "identity", alpha = 0.3) +
    # geom_vline(xintercept = 0) +
    # theme_classic()

