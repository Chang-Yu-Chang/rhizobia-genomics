#' This script plots the trait correlation

library(tidyverse)
library(janitor)
library(cowplot)
library(broom)
# library(lme4) # for linear mixed-effect models
# library(car) # companion to Applied Regression
# #library(factoextra) # for plotting pca ellipses
# library(vegan) # for permanova
source(here::here("analysis/00-metadata.R"))


# Panel A: conceptual plots for the tradeoff concept  ----
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig4A.png"))

# Panel B: one example of permutation ----
boot_pair <- read_csv(paste0(folder_data, "temp/12-01-raw_growth_rate_vs_biomass.csv"))

boot_pair %>%
    #filter(bootstrap %in% 1:40) %>%
    ggplot() +
    geom_point(aes(x = gi0, y = ei0), color = "grey", shape = 21) +
    geom_smooth(aes(x = gi0, y = ei0), color = "grey") +
    geom_point(aes(x = gi1, y = ei1), color = "maroon", shape = 21) +
    geom_smooth(aes(x = gi1, y = ei1), color = "maroon", method = "lm") +
    theme_classic() +
    theme() +
    guides() +
    labs()


boot_pair_across <- select(boot_pair, bootstrap, ends_with("0")) %>% mutate(treatment = "across strain") %>% rename(gi = gi0, ei = ei0)
boot_pair_within <- select(boot_pair, bootstrap, ends_with("1")) %>% mutate(treatment = "within strain") %>% rename(gi = gi1, ei = ei1)

boot_pair_count <- bind_rows(boot_pair_across, boot_pair_within) %>%
    group_by(bootstrap, treatment) %>%
    mutate(across(everything(), rank)) %>%
    ungroup() %>%
    group_by(gi, ei, treatment) %>%
    count()

p2 <- boot_pair_count %>%
    mutate(n = factor(n)) %>%
    ggplot() +
    geom_point(aes(x = gi, y = ei, color = treatment, size = n), shape = 16, alpha = 0.3) +
    scale_color_manual(values = c("across strain" = "steelblue", "within strain" = "maroon")) +
    scale_size_manual(values = setNames(seq(0.1, 3, length.out = 31), 1:31), breaks = c(10,20,30)) +
    facet_grid(~treatment) +
    #scale_x_continuous(expand = c(0,0)) +
    # geom_point(aes(x = gi0, y = ei0), color = "grey", shape = 21) +
    # geom_smooth(aes(x = gi0, y = ei0), color = "grey", method = "lm") +
    # geom_point(aes(x = gi1, y = ei1), color = "maroon", shape = 21) +
    # geom_smooth(aes(x = gi1, y = ei1), color = "maroon", method = "lm") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank()
    ) +
    guides(color = "none") +
    labs(x = "rank[growth rate]", y = "rank[shoot biomass]")



# Panel C ----
p3 <- ggdraw() + draw_image(here::here("plots/cartoons/blank.png"))

# Panel D ----
p4 <- ggdraw() + draw_image(here::here("plots/cartoons/blank.png"))

p <- plot_grid(p1, p2, p3, p4, nrow = 2, axis = "t", align = "h", labels = LETTERS[1:4], scale = 0.9, rel_widths = c(1, 2)) + paint_white_background()

ggsave(here::here("plots/Fig4.png"), p, width = 8, height = 6)













