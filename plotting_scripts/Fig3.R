#' This script

library(tidyverse)
library(cowplot)
library(car) # companion to Applied Regression
library(vegan) # for permanova
library(lme4)
library(emmeans) # estimate marginal means
library(effectsize) # companion to Applied Regression
library(vcd) # for computing effect sizes of categorical response
library(boot)
source(here::here("metadata.R"))

set.seed(1)
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

# PCA plots ----
tbs <- tibble(
    gradient = c("elevation", "urbanization", "elevation", "urbanization", "elevation"),
    exp_plant = c("lupulina", "lupulina", "sativa", "sativa", "sativa"),
    exp_nitrogen = c("N-", "N-", "N-", "N-", "N+")
) %>%
    mutate(treatment = paste(gradient, exp_plant, exp_nitrogen, sep = "\t")) %>%
    mutate(plants_data = list(
        # Lupulinas. First plant exp
        plants1 = plants %>%
            filter(exp_plant == "lupulina", exp_id != "control", gradient == "elevation") %>%
            select(population, exp_id, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
            drop_na(),
        plants2 =  plants %>%
            filter(exp_plant == "lupulina", exp_id != "control", gradient == "urbanization") %>%
            select(population, exp_id, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
            drop_na(),
        # sativas. Second plant exp
        plants3 = plants %>%
            filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", gradient == "elevation") %>%
            select(population, exp_id, shoot_height, nodule_number, longest_petiole_length, leaf_number, leaf_color) %>%
            drop_na(),
        plants4 = plants %>%
            filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", gradient == "urbanization") %>%
            select(population, exp_id, shoot_height, nodule_number, leaf_number, leaf_color, lateral_root_number, longest_lateral_root_length) %>%
            drop_na(),
        plants5 = plants %>%
            filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "with nitrogen", gradient == "elevation") %>%
            select(population, exp_id, shoot_height, nodule_number, longest_petiole_length, leaf_number, leaf_color) %>%
            drop_na()
    ))


do_pca <- function(x) {
    prcomp(x[,-c(1,2)], scale. = TRUE)
}
get_pcs <- function (plants_subset, pca_result) {
    #' Extract the PCs from the pca object
    as_tibble(`[[`(pca_result, "x")) %>%
        mutate(population = plants_subset$population, exp_id = plants_subset$exp_id) %>%
        left_join(distinct(isolates, gradient, population))
}
get_pcvar <- function (pca_result) {
    #' Get the variance explained by the top two PCs
    summary(pca_result)$importance[2, ] %>% round(3) * 100
}
plot_pca <- function (pcsi, pca_result) {
    dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
    # strata by exp_id
    #mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = exp_id))
    mod <- adonis2(dm ~ population, data = pcsi, permutations = 10000)
    pcsi %>%
        ggplot() +
        geom_point(aes(x = PC1, y = PC2, color = population), shape = 21, stroke = 1, size = 2) +
        stat_ellipse(aes(x = PC1, y = PC2, fill = population), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
        annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(pcsi)), hjust = 1.1, vjust = 1.1) +
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
        geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
        geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
        scale_color_manual(values = population_colors) +
        scale_fill_manual(values = population_colors, name = "population") +
        scale_x_continuous(breaks = seq(-8,8,2)) +
        scale_y_continuous(breaks = seq(-8,8,2)) +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey10", fill = NA),
            legend.background = element_blank(),
            plot.background = element_blank(),
            plot.margin = margin(0,0,0,0, "mm")
        ) +
        guides(fill = "none", color = "none") +
        labs(x = paste0("PC1 (", get_pcvar(pca_result)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_result)[2], "%)"))

}

tbs_pca <- tbs %>%
    mutate(
        pca_results = map(plants_data, do_pca),
        pcs = map2(plants_data, pca_results, get_pcs),
        p_pca = map2(pcs, pca_results, plot_pca)
    )
p_pcas <- tbs_pca$p_pca


# Effect size ----
# Remove control and nonsymbiontic strains
plants <- plants %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")
lupulinas <- plants %>% filter(exp_plant == "lupulina")
sativas <- plants %>% filter(exp_plant == "sativa")
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    select(genome_id = genome, species = contig_species, exp_id) %>%
    bind_rows(tibble(genome_id = c("g_src1", "g_bg1"), species = NA, exp_id = c("src-1", "bg-1")))

subset_plants <- function (x, y, z) {
    #' Subset plants data for each trait
    get(paste0(x, "s")) %>% filter(exp_nitrogen == y, gradient == z)
}
compute_cohensd <- function (tb, response) {
    #' Compute Cohen's d
    formu <- paste0(response, " ~ population + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    emm.mod <- emmeans(mod, specs = "population")
    es <- eff_size(emm.mod, sigma = sigma(mod), edf = df.residual(mod))
    return(as_tibble(es))
}

treatments <- tibble(
    pop = c(rep("elevation", 3), rep("urbanization", 3), rep("elevation", 14), rep("urbanization", 7)),
    plant = c(rep("lupulina", 6), rep("sativa", 21)),
    nt = c(rep("without nitrogen", 6+7), rep("with nitrogen", 7), rep("without nitrogen", 7)),
    response = c(rep(c("nodule_number", "shoot_biomass_mg", "root_biomass_mg"), 2),
                 rep(c("nodule_number", "primary_root_nodule_number", "lateral_root_nodule_number", "shoot_height", "longest_petiole_length", "leaf_number", "leaf_color"), 2),
                 "nodule_number", "shoot_height", "leaf_number", "leaf_color", "primary_root_length", "lateral_root_number", "longest_lateral_root_length")
)

# Compute effect size for each trait
treatments_eff <- treatments %>%
    mutate(
        plants_data = pmap(list(plant, nt, pop), subset_plants),
        cohensd = map2(plants_data, response, compute_cohensd)
    )

# Clean the names
ess <- treatments_eff %>%
    unnest(cohensd) %>%
    clean_names() %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "N-",
        nt == "with nitrogen" ~ "N+"
    )) %>%
    mutate(
        pop = factor(pop, c("elevation", "urbanization")),
        plant = factor(plant, c("lupulina", "sativa"))
    ) %>%
    mutate(response = str_remove(response, "mg") %>% str_replace_all("_", " ")) %>%
    arrange(pop, plant, response)

background_df <- tibble(
    pop = factor(c("elevation", "urbanization", "elevation", "urbanization", "elevation"), c("elevation", "urbanization")),
    plant = c("lupulina", "lupulina", "sativa", "sativa", "sativa"),
    host_type = c("source", "source", "alternative", "alternative", "alternative"),
    nt = c("N-", "N-", "N-", "N-", "N+")
)

clean_trait_names <- function (x) {
    str_split(x, pattern = " ")[[1]] %>% str_sub(1,1) %>% paste(collapse = "") %>% toupper() %>% str_pad(width = 4, side = "right", pad = " ")
}
traits <- tibble(
    response = c("nodule number", "root biomass ", "shoot biomass ", "lateral root nodule number", "leaf color", "leaf number", "longest petiole length", "primary root nodule number", "shoot height", "lateral root number", "longest lateral root length", "primary root length"),
    response_abbr = map_chr(response, clean_trait_names)
)

plot_eff <- function (e, pp, pl, nn) {
    bdf <- background_df %>% filter(pop == pp, plant == pl, nt == nn)

    e %>%
        filter(pop == pp, plant == pl, nt == nn) %>%
        left_join(traits) %>%
        ggplot() +
        #geom_rect(data = bdf, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
        #geom_text(data = bdf, aes(label = paste("M.", plant)), x = Inf, y = -Inf, size = 5, hjust = 1, vjust = -0.4, fontface = "italic") +
        geom_text(data = bdf, aes(label = paste(host_type, "host", nt)), x = Inf, y = -Inf, size = 4, hjust = 0.5, vjust = -0.4) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(x = response_abbr, y = effect_size, shape = nt), size = 3, stroke = 1, fill = "white") +
        geom_linerange(aes(x = response_abbr, ymin = lower_cl, ymax = upper_cl), color = "grey10", linewidth = 1) +
        geom_point(aes(x = response_abbr, y = effect_size, shape = nt), size = 3, stroke = 1, fill = "white") +
        scale_x_discrete(expand = c(0,.8), position = "top") +
        scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
        scale_shape_manual(values = c(`with nitrogen` = 16, `without nitrogen` = 21), labels = c("N+", "N-")) +
        scale_fill_manual(values = plant_colors) +
        #facet_grid(plant~., space = "free_y", switch = "y") +
        coord_flip(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
            panel.grid.minor.y = element_blank(),
            axis.title.y = element_blank(),
            #axis.text.y = element_blank(),
            strip.text = element_blank(),
            legend.title = element_blank(),
            legend.position = "inside",
            legend.position.inside = c(0.13,0.45),
            legend.margin = margin(0,0,0,0),
            legend.key = element_blank(),
            legend.key.spacing = unit(-1, "mm"),
            legend.background = element_rect(color = "grey90", fill = "white"),
            plot.margin = margin(0,10,0,0, "mm"),
            plot.background = element_rect(color = NA, fill = "white")
        ) +
        guides(color = "none", fill = "none", shape = "none") +
        labs(y = "standardized mean difference")
}
p_effs <- list(
    plot_eff(ess, "elevation", "lupulina", "N-") + theme(axis.title.x = element_blank()),
    plot_eff(ess, "urbanization", "lupulina", "N-") + theme(axis.title.x = element_blank()),
    plot_eff(ess, "elevation", "sativa", "N-") + theme(axis.title.x = element_blank()),
    plot_eff(ess, "urbanization", "sativa", "N-"),
    plot_eff(ess, "elevation", "sativa", "N+")
)

# Combine figures ----
p_dat <- plot_grid(
    p_pcas[[1]], p_effs[[1]],
    p_pcas[[2]], p_effs[[2]],
    p_pcas[[3]], p_effs[[3]],
    p_pcas[[4]], p_effs[[4]],
    p_pcas[[5]], p_effs[[5]],
    rel_widths = c(1,2), nrow = 5, scale = .85,
    #rel_heights = c(.1,1,1,1),
    labels = c("A", "", "B", "", "C", "", "D", "", "E", "")
)


p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig3.png"), scale = 1) +
    draw_plot(p_dat, x = .15, y = 0.0, width = .85, height = .96) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

#ggsave(here::here("plots/Fig3.png"), p, width = 8, height = 12)

























