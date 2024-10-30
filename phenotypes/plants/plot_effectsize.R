#' Plot effect size

library(tidyverse)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
cohensds <- read_csv(paste0(folder_data, "phenotypes/effectsize/cohensds.csv"))

# 1. Cohen's d ----
plot_cohensds <- function (cohensds, plant) {
    cohensds %>%
        filter(exp_plant == plant) %>%
        mutate(exp_nitrogen = factor(exp_nitrogen, c("without nitrogen", "with nitrogen"))) %>%
        ggplot() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_linerange(aes(x = trait, ymin = lower_cl, ymax = upper_cl), color = "grey10", linewidth = 1, position = position_dodge2(width = .5)) +
        geom_point(aes(x = trait, y = effect_size, shape = exp_nitrogen), size = 3, stroke = 1, fill = "white", position = position_dodge2(width = .5, reverse = T)) +
        scale_x_discrete(expand = c(0,.8), position = "top") +
        scale_y_continuous(limits = c(-3, 3), expand = c(0,.1), breaks = -3:3) +
        scale_shape_manual(values = c(`with nitrogen` = 16, `without nitrogen` = 21), labels = c("N+", "N-")) +
        scale_fill_manual(values = plant_colors) +
        facet_grid(gradient~exp_nitrogen, scales = "free_y", space = "free_y", switch = "y") +
        coord_flip(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            plot.background = element_rect(color = NA, fill = "white"),
        ) +
        guides(color = "none", fill = "none") +
        labs(y = "standardized mean difference (Cohen's d)", title = paste0(plant, " experiment"))
}
p <- plot_cohensds(cohensds, "sativa")
ggsave(paste0(folder_phenotypes, "effectsize/01-cohensd_sativa.png"), p, width = 8, height = 5)
p <- plot_cohensds(cohensds, "lupulina")
ggsave(paste0(folder_phenotypes, "effectsize/02-cohensd_lupulina.png"), p, width = 5, height = 5)

# 2. Plot the boxplot of traits with Cohen's d -----
set.seed(1)
plot_boxes <- function (plants, gra, plant) {
    strips = strip_nested(
        # Vertical strips
        background_y = elem_list_rect(
            color = NA,
            fill = c("white")
        ),
        text_y = elem_list_text(
            angle = 0
        ),
        bleed = T,
        by_layer_y = F,
        clip = "off", size = "variable"
    )
    plants %>%
        filter(gradient == gra) %>%
        filter(population != "control", exp_plant == plant, exp_nitrogen == "without nitrogen") %>%
        select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
        group_by(gradient, population, exp_plant) %>%
        filter(nodule_number <100) %>%
        pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
        left_join(traits) %>%
        arrange(trait_type) %>%
        ggplot(aes(x = population, y = value, fill = population)) +
        geom_boxplot(outlier.shape = -1, alpha = 0.5) +
        geom_jitter(size = .5, shape = 21, width = .2) +
        scale_fill_manual(values = population_colors) +
        coord_flip(clip = "off") +
        #facet_nested_wrap(~trait_type + trait, scales = "free", ncol = 1, nest_line = element_line(color = "blue", linetype = 2), dir = "v", strip.position = "left", strip = strips) +
        facet_nested_wrap(trait_type + trait_pre ~., ncol = 1, strip.position = "left", axes = "all", scales = "free", solo_line = T, nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1), strip = strips) +
        #facet_wrap2(trait_pre~., scales = "free", remove_labels = F, axes = "all", ncol = 1) +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(color = NA, fill = "grey95"),
            strip.background = element_rect(color = NA),
            #strip.text = element_blank(),
            legend.position = "top",
            legend.key = element_rect(fill = NA, color = NA),
            legend.key.height = unit(10, "mm"),
            legend.text = element_text(size = 8),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.box.margin = unit(c(0,0,-5,0), "mm")
        ) +
        guides(fill = guide_legend(override.aes = list(color = NA, size = 0, shape = 0))) +

         labs()
}
plot_eff <- function (cohends, gra, plant) {
    cohensds %>%
        filter(gradient == gra) %>%
        filter(exp_plant == plant, exp_nitrogen == "without nitrogen") %>%
        mutate(exp_nitrogen = factor(exp_nitrogen, c("without nitrogen", "with nitrogen"))) %>%
        left_join(traits) %>%
        ggplot() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_linerange(aes(x = trait_pre, ymin = lower_cl, ymax = upper_cl), color = "grey10", linewidth = 1, position = position_dodge2(width = .5)) +
        geom_point(aes(x = trait_pre, y = effect_size, shape = exp_nitrogen), size = 3, stroke = 1, fill = "white", position = position_dodge2(width = .5, reverse = T)) +
        scale_x_discrete(expand = c(0,.8), position = "top") +
        scale_y_continuous(limits = c(-3, 3), expand = c(0,.1), breaks = -3:3) +
        scale_shape_manual(values = c(`with nitrogen` = 16, `without nitrogen` = 21), labels = c("N+", "N-")) +
        scale_fill_manual(values = plant_colors) +
        #facet_nested_wrap(gradient + trait~., remove_labels = "x", axes = "all", scales = "free_y", dir = "v", strip.position = "left", ncol = 1, strip = strips) +
        #facet_wrap2(vars(trait_pre), scales = "free_y", ncol = 1, axes = "all", remove_labels = "x") +
        #facet_nested(trait_type + trait_pre ~., scales = "free_y", axes = "all", remove_labels = "x", nest_line = element_line(color = "grey30", linetype = 2), strip = strips) +
        facet_wrap2(trait_pre~., scales = "free_y", remove_labels = "x", axes = "all", ncol = 1) +
        coord_flip(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(color = NA, fill = "grey95"),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            strip.background = element_rect(color = NA),
            strip.text = element_blank(),
            #strip.text = element_text(angle = 90),
            plot.background = element_rect(color = NA, fill = "white")
        ) +
        guides(color = "none", fill = "none") +
        labs(y = "standardized mean difference")
}

p1 <- plot_boxes(plants, "elevation", "sativa")
p2 <- plot_eff(plants, "elevation", "sativa")
p <- plot_grid(p1, p2, align = "h", axis = "tb", ncol = 2, rel_widths = c(1, .5))
ggsave(paste0(folder_phenotypes, "effectsize/03-elevation.png"), p, width = 8, height = 5)

p1 <- plot_boxes(plants, "urbanization", "sativa")
p2 <- plot_eff(plants, "urbanization", "sativa")
p <- plot_grid(p1, p2, align = "h", axis = "tb", ncol = 2, rel_widths = c(1, .5))
ggsave(paste0(folder_phenotypes, "effectsize/04-urbanization.png"), p, width = 8, height = 5)

# Combined ----
p1 <- plot_boxes(plants, "elevation", "sativa")
p2 <- plot_eff(plants, "elevation", "sativa")
p3 <- plot_boxes(plants, "urbanization", "sativa")
p4 <- plot_eff(plants, "urbanization", "sativa")
p <- plot_grid(p1, p2, p3, p4, align = "h", axis = "tb", ncol = 2, rel_widths = c(1, .5), labels = LETTERS[1:4]) +
    draw_text(c("Elevation", "Urbanization"), x = 0.05, y = c(0.98,  0.48), hjust = 0)
    #draw_text(c("trait type", "trait"), x = c(0.05, 0.16), y = 0.45, hjust = 0.5, size = 10)
ggsave(paste0(folder_phenotypes, "effectsize/05-traits_effectsize.png"), p, width = 8, height = 10)




