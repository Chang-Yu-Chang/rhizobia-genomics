#' This script plots the pca figs

library(tidyverse)
library(cowplot)
library(vegan) # for permanova
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

# Master table for subset
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
do_permanova <- function (pcs) {
    #dm <- vegdist(select(pcs, starts_with("PC")), method = "euclidean")
    dm <- dist(select(pcs, starts_with("PC")))
    # strata by populaiton
    mod <- adonis2(dm ~ population, data = pcs, permutations = 1000)
    return(mod)
}
plot_pca <- function (pcsi, pca_result, permanova_result) {
    signlab <- clean_p_lab(permanova_result[1,5])
    ann <- paste0("N_rhi=", length(unique(pcsi$exp_id)),", N_plant=", nrow(pcsi), ", permanova: ", signlab)

    pcsi %>%
        #left_join(select(isolates, exp_id, genome_id)) %>%
        ggplot() +
        stat_ellipse(aes(x = PC1, y = PC2, group = exp_id, fill = population), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
        geom_point(aes(x = PC1, y = PC2, color = population), shape = 21, stroke = 1, size = 2) +
        # annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(pcsi)), hjust = 1.1, vjust = 1.1) +
        # annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.6, label = signlab) +
        #annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.6, label = ann) +
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
            legend.position = "none"
        ) +
        guides(color = "none") +
        labs(subtitle = ann, x = paste0("PC1 (", get_pcvar(pca_result)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_result)[2], "%)"))

}

pcsi <- tbs_pca$pcs[[3]]
pca_result <- tbs_pca$pca_result[[3]]
permanova_result <- tbs_pca$permanova_result[[3]]
plants_data <- tbs_pca$plants_data[[3]]


tbs_pca <- tbs %>%
    mutate(
        pca_result = map(plants_data, do_pca),
        pcs = map2(plants_data, pca_result, get_pcs),
        permanova_result = map(pcs, do_permanova),
        p_pca = pmap(list(pcs, pca_result, permanova_result), plot_pca)
    )

p <- plot_grid(
    tbs_pca$p_pca[[1]] + ggtitle("elevation lupulina N-"),
    tbs_pca$p_pca[[2]] + ggtitle("urbanization lupulina N-"),
    tbs_pca$p_pca[[3]] + ggtitle("elevation sativa N-"),
    tbs_pca$p_pca[[4]] + ggtitle("urbanization sativa N-"),
    tbs_pca$p_pca[[5]] + ggtitle("elevation sativa N+"),
    ncol = 2, scale = .95, align = "vh", axis = "lr") +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_phenotypes, "pca/01-pca.png"), p, width = 9, height = 12)


