#' This script compares the symbiosis traits and effect size
#' 1. Prepare data
#' 2. Check assumptions
#' 3. Run models
#' 4. Make the stat tables
#' 5. Plot

library(tidyverse)
library(cowplot)
library(ggh4x)
library(grid)
library(vegan) # for permanova
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
cohensds <- read_csv(paste0(folder_data, "phenotypes/plants/effectsize/cohensds.csv"))

# 1. Prepare data ----
plants_n <- plants %>%
    filter(population != "control", exp_nitrogen == "N-") %>%
    select(
        -nodule_shape, -nodule_size, -nodule_color, -exp_labgroup, -exp_labsection,
        -primary_root_nodule_number, -lateral_root_nodule_number,
        -longest_petiole_length, -longest_lateral_root_length,
        -lateral_root_number, -primary_root_length
    ) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    left_join(isolates) %>%
    ungroup()

cohensds_n <- cohensds %>%
    filter(!str_detect(trait, "primary|lateral|root_length|longest"))

# 2. Plots ----
plot_boxes <- function (plants_n, gra, plant, nt) {
    # gra = "elevation"
    # plant = "sativa"
    # gra = "urbanization"

    # Order the strains
    exp_id_lev <- plants_n %>% distinct(exp_id, .keep_all = T) %>%
        mutate(population = factor(population, rev(c("high elevation", "low elevation", "suburban", "urban")))) %>%
        arrange(population) %>% pull(exp_id)

    plants_n %>%
        filter(gradient == gra) %>%
        filter(exp_plant == plant) %>%
        filter(exp_nitrogen == nt) %>%
        left_join(traits) %>%
        mutate(trait_type = factor(trait_type, c("shoot", "nodule", "leaf", "root"))) %>%
        mutate(trait_pre2 = factor(trait_pre2, traits$trait_pre2)) %>%
        mutate(exp_id = factor(exp_id, exp_id_lev)) %>%
        arrange(trait_type, trait_pre2) %>%
        group_by(gradient, population, exp_plant, trait_type, trait_pre2, value, exp_id) %>%
        count() %>%
        ggplot() +
        geom_boxplot(aes(x = exp_id, y = value, fill = population), alpha = 0.3, width = .7, outlier.size = -1) +
        geom_point(aes(x = exp_id, y = value, color = population, size = n), alpha = .4, shape = 16, position = position_dodge(width = .7)) +
        scale_fill_manual(values = population_colors) +
        scale_color_manual(values = population_colors) +
        scale_size_continuous(range = c(.5,3), limits = c(1, 20), breaks = c(5, 10 ,20)) +
        coord_flip(clip = "off") +
        facet_nested_wrap(
            trait_type + trait_pre2 ~., ncol = 1, strip.position = "left", axes = "all", scales = "free",
            solo_line = T, nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1),
            strip = strip_nested(
                background_y = elem_list_rect(color = NA, fill = c("white")),
                text_y = elem_list_text(size = 8,angle = 0),
                bleed = T, by_layer_y = F, clip = "off", size = "variable"
            )) +
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
            legend.text = element_text(size = 6),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.box.margin = unit(c(0,0,-5,0), "mm"),
            plot.background = element_blank()
        ) +
        guides(fill = guide_legend(override.aes = list(color = NA, size = 0, shape = 0)), color = "none") +
        labs()
}
plot_cohensds <- function (cohensds, gra, plant, nt) {
    cohensds %>%
        filter(gradient == gra) %>%
        filter(exp_plant == plant) %>%
        filter(exp_nitrogen == nt) %>%
        mutate(exp_nitrogen = factor(exp_nitrogen, c("N-", "N+"))) %>%
        left_join(traits) %>%
        mutate(trait_pre2 = factor(trait_pre2, rev(levels(trait_pre2)))) %>%
        arrange(trait_type) %>%
        group_by(gradient, exp_plant, trait_type, trait_pre2) %>%
        ggplot() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_linerange(aes(x = trait_pre2, ymin = lower_cl, ymax = upper_cl), color = "grey10", linewidth = 1, position = position_dodge2(width = .5)) +
        geom_point(aes(x = trait_pre2, y = effect_size, shape = exp_nitrogen), size = 3, stroke = 1, fill = "white", position = position_dodge2(width = .5, reverse = T)) +
        #facet_grid2(trait_pre2~., scales = "free_y") +
        facet_wrap2(trait_pre2 ~., ncol = 1, strip.position = "left", axes = "all", scales = "free_y", remove_labels = "all") +
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
p7 <- plot_boxes(plants_n, "urbanization", "lupulina", "N-")
leg1 <- get_legend(p1 + theme(legend.position = "right", legend.spacing = unit(-8, "mm"), legend.key.size = unit(6, "mm"), legend.text = element_text(size = 10)) + guides(size = "none", fill = guide_legend(order = 2, nrow = 1, override.aes = list(color = NA))))
leg2 <- get_legend(p1 + theme(legend.position = "right", legend.spacing = unit(-8, "mm"), legend.key.size = unit(6, "mm"), legend.text = element_text(size = 10)) + guides(size = guide_legend(nrow = 1, order = 1), fill = "none"))
leg3 <- get_legend(p3 + theme(legend.position = "right", legend.spacing = unit(-8, "mm"), legend.key.size = unit(6, "mm"), legend.text = element_text(size = 10)) + guides(size = "none", fill = guide_legend(order = 2, nrow = 1, override.aes = list(color = NA))))

p2 <- plot_cohensds(cohensds_n, "elevation", "sativa", "N-") + theme(axis.title.x = element_blank())
p4 <- plot_cohensds(cohensds_n, "urbanization", "sativa", "N-") + theme(axis.title.x = element_blank())
p6 <- plot_cohensds(cohensds_n, "elevation", "lupulina", "N-")
p8 <- plot_cohensds(cohensds_n, "urbanization", "lupulina", "N-")


arrow_y = 0.92
arrow_grob1 <- linesGrob(x = unit(c(.35, .48),"npc"), y = unit(c(arrow_y, arrow_y), "npc"), gp = gpar(col = "black", lwd = 2), arrow =  arrow(length = unit(2, "mm"), angle = 30, type = "open", ends = "both"))
arrow_grob2 <- linesGrob(x = unit(c(.85, .98),"npc"), y = unit(c(arrow_y, arrow_y), "npc"), gp = gpar(col = "black", lwd = 2), arrow =  arrow(length = unit(2, "mm"), angle = 30, type = "open", ends = "both"))

p_combined <- plot_grid(
    leg1, leg2, leg3, NULL,
    p1 + theme(legend.position = "none"), p2, p3 + theme(legend.position = "none"), p4,
    p5 + theme(legend.position = "none"), p6, p7 + theme(legend.position = "none"), p8,
    ncol = 4,  align = "vh", axis = "tbr", scale = .95,
    rel_widths = c(1,.5,1,.5), rel_heights = c(.5,4,3),
    labels = c(rep("", 4), "A", "", "B", "", rep("", 0), "C", "", "D", ""), label_x = .1
) +
    # Panel titles legend
    draw_text("N strains = 8, N plants = 225", x = .1, y = .91, size = 10, hjust = 0, vjust = -1) +
    draw_text("N strains = 8, N plants = 105", x = .6, y = .91, size = 10, hjust = 0, vjust = -1) +
    draw_text("N strains = 6, N plants = 159", x = .1, y = .38, size = 10, hjust = 0, vjust = -1) +
    draw_text("N strains = 8, N plants = 68", x = .6, y = .38, size = 10, hjust = 0, vjust = -1) +
    # Effect size legend
    draw_grob(arrow_grob1) +
    draw_grob(arrow_grob2) +
    draw_text("low", x = .34, y = arrow_y, size = 10, hjust = 0, vjust = -1) +
    draw_text("high", x = .49, y = arrow_y, size = 10, hjust = 1, vjust = -1) +
    draw_text("urban", x = .84, y = arrow_y, size = 10, hjust = 0, vjust = -1) +
    draw_text("suburban", x = .99, y = arrow_y, size = 10, hjust = 1, vjust = -1)

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig2.png"), scale = 1) +
    draw_plot(p_combined, x = .05, width = .9) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/Fig2.png"), p, width = 10, height = 10)


#
plants %>% filter(exp_nitrogen == "N-", exp_id != "control") %>% group_by(exp_plant, gradient) %>% count()

# 3. PERMANOVA ----
set.seed(1)
# Sativa elevation
dat <- plants %>%
    filter(population != "control", exp_nitrogen == "N-") %>%
    filter(gradient == "elevation", exp_plant == "sativa") %>%
    drop_na(shoot_height, nodule_number, leaf_color, leaf_number) %>%
    select(gradient, population, site, exp_id, exp_waterblock, shoot_height, nodule_number, leaf_color, leaf_number)
m <- select(dat, shoot_height, nodule_number, leaf_color, leaf_number)
adonis2(m ~ population, data = dat, permutation = 1000)

# Sativa urbanization
dat <- plants %>%
    filter(population != "control", exp_nitrogen == "N-") %>%
    filter(gradient == "urbanization", exp_plant == "sativa") %>%
    drop_na(shoot_height, nodule_number, leaf_color, leaf_number) %>%
    select(gradient, population, site, exp_id, exp_waterblock, shoot_height, nodule_number, leaf_color, leaf_number)
m <- select(dat, shoot_height, nodule_number, leaf_color, leaf_number)
adonis2(m ~ population, data = dat, permutation = 1000)

# Luplina elevation
dat <- plants %>%
    filter(population != "control", exp_nitrogen == "N-") %>%
    filter(gradient == "elevation", exp_plant == "lupulina") %>%
    drop_na(shoot_biomass_mg, nodule_number, root_biomass_mg) %>%
    select(gradient, population, site, exp_id, exp_waterblock, shoot_biomass_mg, nodule_number, root_biomass_mg)
m <- select(dat, shoot_biomass_mg, nodule_number, root_biomass_mg)
adonis2(m ~ population, data = dat, permutation = 1000)

# Luplina urbnaization
dat <- plants %>%
    filter(population != "control", exp_nitrogen == "N-") %>%
    filter(gradient == "urbanization", exp_plant == "lupulina") %>%
    drop_na(shoot_biomass_mg, nodule_number, root_biomass_mg) %>%
    select(gradient, population, site, exp_id, exp_waterblock, shoot_biomass_mg, nodule_number, root_biomass_mg)
m <- select(dat, shoot_biomass_mg, nodule_number, root_biomass_mg)
adonis2(m ~ population, data = dat, permutation = 1000)
