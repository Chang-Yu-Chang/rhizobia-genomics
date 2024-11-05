#' This script compares the symbiosis traits and effect size
#' 1. Prepare data
#' 2. Check assumptions
#' 3. Run models
#' 4. Make the stat tables
#' 5. Plot

library(tidyverse)
library(cowplot)
library(ggh4x)
# library(flextable)
# library(broom.mixed) # for tidying the model outputs
# library(lme4) # for lmer
# library(car) # for anova
# library(boot) # for bootstrapping
source(here::here("metadata.R"))
#options(contrasts=c("contr.sum", "contr.poly"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
cohensds <- read_csv(paste0(folder_data, "phenotypes/effectsize/cohensds.csv"))

# 1. Prepare data ----
plants_n <- plants %>%
    filter(population != "control", exp_nitrogen == "N-") %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
    select(-primary_root_nodule_number, -lateral_root_nodule_number) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    left_join(isolates) %>%
    ungroup()

cohensds_n <- cohensds %>%
    filter(!trait %in% c("primary_root_nodule_number", "lateral_root_nodule_number"))

# 2. Plots ----
plot_boxes <- function (plants_n, gra, plant, nt) {
    # gra = "elevation"
    # plant = "sativa"
    strips = strip_nested(
        background_y = elem_list_rect(
            color = NA,
            fill = c("white")
        ),
        text_y = elem_list_text(
            size = 8,
            angle = 0
        ),
        bleed = T,
        by_layer_y = F,
        clip = "off", size = "variable"
    )
    # ttp <- tb_tidied_p %>%
    #     filter(gradient == gra)


    plants_n %>%
        filter(gradient == gra) %>%
        filter(exp_plant == plant) %>%
        filter(exp_nitrogen == nt) %>%
        left_join(traits) %>%
        mutate(population = factor(population, c("low elevation", "high elevation", "urban", "suburban"))) %>%
        arrange(trait_type) %>%
        group_by(gradient, population, exp_plant, trait_type, trait_pre, value, exp_id) %>%
        count() %>%
        ggplot() +
        geom_boxplot(aes(x = population, y = value, fill = population), alpha = 0.3, width = .7, outlier.size = -1) +
        geom_point(aes(x = population, group = exp_id, y = value, color = population, size = n), alpha = .4, shape = 16, position = position_dodge(width = .7)) +
        # stat_summary(
        #     fun.data = function(x) {
        #         y_min <- min(x)
        #         y_max <- max(x)
        #         y_segment <- y_min + 1.05 * (y_max - y_min)  # 75% line
        #         data.frame(y = y_segment)
        #     },
        #     geom = "segment",
        #     aes(x = 1, xend = 2, y = value, yend = value),
        #     color = "black", linetype = 1, linewidth = .5
        # ) +
        #geom_text(data = ttp, aes(label = signlab), x = 1.5, y = Inf, hjust = 0.3) +
        scale_fill_manual(values = population_colors) +
        scale_color_manual(values = population_colors) +
        scale_size_continuous(range = c(.5,3)) +
        coord_flip(clip = "off") +
        facet_nested_wrap(trait_type + trait_pre ~., ncol = 1, strip.position = "left", axes = "all", scales = "free", solo_line = T, nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1), strip = strips) +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(color = NA, fill = "grey95"),
            strip.background = element_rect(color = NA),
            legend.position = "top",
            legend.key = element_rect(fill = NA, color = NA),
            legend.key.height = unit(10, "mm"),
            legend.text = element_text(size = 8),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.box.margin = unit(c(0,0,-5,0), "mm")
        ) +
        guides(size = "none", fill = guide_legend(override.aes = list(color = NA, size = 0, shape = 0), reverse = T), color = "none") +
        labs()
}
plot_cohensds <- function (cohensds, gra, plant, nt) {
    cohensds %>%
        filter(gradient == gra) %>%
        filter(exp_plant == plant) %>%
        filter(exp_nitrogen == nt) %>%
        mutate(exp_nitrogen = factor(exp_nitrogen, c("N-", "N+"))) %>%
        left_join(traits) %>%
        mutate(trait_pre = factor(trait_pre, rev(levels(trait_pre)))) %>%
        arrange(trait_type) %>%
        group_by(gradient, exp_plant, trait_type, trait_pre) %>%
        ggplot() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_linerange(aes(x = trait_pre, ymin = lower_cl, ymax = upper_cl), color = "grey10", linewidth = 1, position = position_dodge2(width = .5)) +
        geom_point(aes(x = trait_pre, y = effect_size, shape = exp_nitrogen), size = 3, stroke = 1, fill = "white", position = position_dodge2(width = .5, reverse = T)) +
        #facet_grid2(trait_pre~., scales = "free_y") +
        facet_wrap2(trait_pre ~., ncol = 1, strip.position = "left", axes = "all", scales = "free_y", remove_labels = "all") +
        scale_x_discrete(expand = c(0,.8), position = "top") +
        scale_y_continuous(limits = c(-3, 3), expand = c(0,0), breaks = -3:3) +
        scale_fill_manual(values = plant_colors) +
        coord_flip(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            #panel.grid.major.y = element_blank(),
            panel.background = element_rect(color = NA, fill = "grey95"),
            panel.border = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            strip.background = element_rect(color = NA),
            strip.text = element_blank(),
            plot.background = element_rect(color = NA, fill = "white"),
        ) +
        guides(color = "none", fill = "none") +
        labs(y = "standardized mean difference")
}

p1 <- plot_boxes(plants_n, "elevation", "sativa", "N-")
p3 <- plot_boxes(plants_n, "urbanization", "sativa", "N-")
p5 <- plot_boxes(plants_n, "elevation", "lupulina", "N-")
leg1 <- get_legend(p1 + theme(legend.position = "right", legend.key.size = unit(1, "cm"), legend.text = element_text(size = 10)))
leg3 <- get_legend(p3 + theme(legend.position = "right", legend.key.size = unit(1, "cm"), legend.text = element_text(size = 10)))


p2 <- plot_cohensds(cohensds_n, "elevation", "sativa", "N-")
p4 <- plot_cohensds(cohensds_n, "urbanization", "sativa", "N-")
p6 <- plot_cohensds(cohensds_n, "elevation", "lupulina", "N-")


p_left <- plot_grid(
    p1 + theme(legend.position = "none"),
    p2,
    p5 + theme(legend.position = "none"),
    p6,
    ncol = 2, labels = c("A", "", "C", ""), scale = 0.98,
    align = "vh", axis = "tbr", rel_widths = c(1, .5), rel_heights = c(1, .5)
)
p_right <- plot_grid(
    p3 + theme(legend.position = "none"),
    p4, leg1, leg3,
    ncol = 2, labels = c("B", "", "", ""), scale = 0.98,
    align = "h", axis = "tb", rel_widths = c(1, .5), rel_heights = c(1, .2)
)

p <- plot_grid(p_left, p_right, nrow = 1, align = "vh", axis = "tbr") +
    theme(plot.background = element_rect(fill = "white", color = NA))


# p <- plot_grid(
#     p1, p2, p3, p4,
#     p5 + theme(legend.position = "none"),
#     p6, NULL, NULL,
#     ncol = 4, labels = c("A", "", "B", "", "C", ""), scale = 0.98,
#     align = "vh", axis = "tbr", rel_widths = c(1, .5, 1, .5), rel_heights = c(1, .5, 1, .5)
# ) +
#     theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_phenotypes, "plants/symbiosis.png"), p, width = 10, height = 8)

