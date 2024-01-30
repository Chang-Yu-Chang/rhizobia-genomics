#' This script plots the trait data of plant inoculation experiments

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
library(ellipse) # for calculating the ellipse
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
sysplants <- read_csv(paste0(folder_data, "temp/23-plants.csv"))
plants_long <- read_csv(paste0(folder_data, "temp/23-plants_long.csv"))
plants_wide <- read_csv(paste0(folder_data, "temp/23-plants_wide.csv"))


# 1. plot all plant traits
length(unique(plants_long$id)) # 231 plants
length(unique(plants_long$exp_id)) # 14 rhizobia + 1 control
plants_long %>%
    arrange(site_group) %>%
    distinct(id, site_group, exp_id) %>%
    group_by(site_group, exp_id) %>%
    count()

p <- plants_long %>%
    ggplot() +
    geom_boxplot(aes(x = exp_id, y = value, color = site_group), outlier.shape = NA) +
    geom_jitter(aes(x = exp_id, y = value, color = site_group), width = 0.2, height = 0, shape = 21) +
    facet_grid(trait~site_group, scales = "free", space = "free_x") +
    #scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    guides(color = "none") +
    labs()

ggsave(paste0(folder_data, "temp/23a-01-plant_traits.png"), p, width = 10, height = 10)

# 2. Plot the nodule vs. plant biomass
p <- plants_wide %>%
    filter(site_group != "control") %>%
    ggplot() +
    geom_point(aes(x = nodule_number, y = dry_weight_mg, color = site_group), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
    scale_color_manual(values = brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)]) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    coord_equal() +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey95")
    ) +
    guides() +
    labs(x = "# of nodules", y = "shoot biomass (mg)")

ggsave(paste0(folder_data, "temp/23a-02-nodule_vs_shoot.png"), p, width = 5, height = 4)


# 3. facet plot ----
p <- plants_wide %>%
    filter(site_group != "control") %>%
    ggplot() +
    geom_point(aes(x = nodule_number, y = dry_weight_mg, color = site_group), shape = 21, size = 2, stroke = 1)+
    scale_color_manual(values = brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)]) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    facet_wrap(~site_group, nrow = 2) +
    coord_equal() +
    theme_classic() +
    #theme_facets +
    # theme(
    #     panel.grid.major = element_line(color = "grey95"),
    #     strip.background = element_rect(fill = "grey95", color = NA),
    #     panel.border = element_rect(fill = NA, color = "black")
    # ) +
    guides() +
    labs(x = "# of nodules", y = "shoot biomass (mg)")

ggsave(paste0(folder_data, "temp/23a-03-nodule_vs_shoot_facet.png"), p, width = 8, height = 6)


"Work on making the correlation within composite traits"


















# Stat for slope in regression
lm_x_y <- function (dat, x_var, y_var) {
    formula <- as.formula(paste(y_var, "~", x_var))
    lm(formula, data = dat) %>%
        broom::tidy() %>%
        filter(term == x_var) %>%
        select(term, estimate, p.value)
}

plants_wide %>%
    filter(site_group != "control") %>%
    nest(data = c(id, exp_id, dry_weight_mg, nodule_number, root_weight_mg)) %>%
    mutate(lm = map(data, ~lm_x_y(.x, "nodule_number", "dry_weight_mg"))) %>%
    select(site_group, lm) %>%
    unnest(lm)

# 4. plot the root traits ----
p <- plants %>%
    #select(id, exp_id, site_group, all_of(traits2)) %>%
    #select(id, site_group, exp_id, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>%
    pivot_longer(cols = c(-id, -exp_id, -site_group), names_to = "trait") %>%
    ggplot() +
    geom_boxplot(aes(x = exp_id, y = value, color = site_group), outlier.shape = NA) +
    geom_jitter(aes(x = exp_id, y = value, color = site_group), width = 0.2, height = 0, shape = 21) +
    facet_grid(trait~site_group, scales = "free", space = "free_x") +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme_facets +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    guides(color = "none") +
    labs()

ggsave(paste0(folder_data, "temp/23a-04-root_traits.png"), p, width = 8, height = 40)

# 5. plot the trait correlation ----
# Names
plants_traits <- plants %>%
    filter(site_group %in% c("high-elevation", "low-elevation")) %>%
    select(all_of(traits))
tb_traits <- tibble(trait = ordered(names(plants_traits), names(plants_traits))) %>%
    expand(trait, trait) %>%
    rename(trait1 = trait...1, trait2 = trait...2) %>%
    filter(trait1 < trait2) %>%
    arrange(trait1, trait2)

# Calculate cor
cor_matrix <- cor(plants_traits, use = "complete.obs")
cor_x_y <- function (dat) cor.test(unlist(dat[,1]), unlist(dat[,2]), method = "spearman", exact = F) %>% broom::tidy() %>% clean_names()
tb_cortest <- tb_traits %>%
    rowwise() %>%
    mutate(traits = list(select(plants_traits, trait1, trait2))) %>%
    ungroup() %>%
    # Compute cor test
    mutate(cortest = map(traits, ~cor_x_y(.x))) %>%
    select(-traits) %>%
    unnest(cortest) %>%
    select(trait1, trait2, cor=estimate, p_value, method, alternative) %>%
    mutate(asterisk = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        T ~ "n.s."
    ))

# tb_cor <- as_tibble(cor_matrix) %>%
#     mutate(trait1 = colnames(.)) %>%
#     pivot_longer(-trait1, names_to = "trait2", values_to = "cor") %>%
#     mutate(trait1 = ordered(trait1, names(plants_traits))) %>%
#     mutate(trait2 = ordered(trait2, names(plants_traits))) %>%
#     filter(trait1 < trait2) %>%
#     arrange(trait1, trait2)




# Calculate ellipse
tb_ell <- tb_traits %>%
    rowwise() %>%
    # Calculate ellipse
    mutate(cor_m = list(ellipse(cor_matrix, which = which(names(plants_traits) %in% c(trait1, trait2))))) %>%
    mutate(cor_tb = list(tibble(x = cor_m[,which(colnames(cor_m) == trait1)], y = cor_m[,which(colnames(cor_m) == trait2)]))) %>%
    select(-cor_m) %>%
    unnest(cor_tb)

tb_cortest_swap <- tibble(trait1 = tb_cortest$trait2, trait2 = tb_cortest$trait1, cor = tb_cortest$cor, asterisk = tb_cortest$asterisk)

#
p <- tb_ell %>%
    left_join(tb_cortest) %>%
    ggplot() +
    geom_polygon(aes(x = x, y = y, fill = cor), color = "black") +
    geom_text(data = tb_cortest_swap, x = 0, y = 0, aes(label = round(cor, 2))) +
    geom_text(data = tb_cortest_swap, x = 0, y = 0, aes(label = asterisk), vjust = -0.5) +
    facet_grid(trait2 ~ trait1, switch = "y", drop = F) +
    scale_fill_gradient2(low = "maroon", high = "steelblue", breaks = c(1,0,-1), limits = c(-1,1)) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        #panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_rect(color = "white", fill = "white")
    ) +
    guides() +
    labs()


ggsave(paste0(folder_data, "temp/23a-05-trait_cor.png"), p, width = 12, height = 10)












