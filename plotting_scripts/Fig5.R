#' This script plots the tree

renv::load()
library(tidyverse)
library(cowplot)
library(tidytree)
library(ggtree)
#library(ggrepel)
library(ggh4x) # for nested strips
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates_traits <- read_csv(paste0(folder_data, "phenotypes_analysis/isolates_traits.csv"))
traits_ps <- read_csv(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))
traits_ps <- traits_ps %>% mutate(across(c(k, p_k, lambda, p_lambda), function(x) round(x, 3)))
isolates <- isolates %>%
    left_join(isolates_contigs) %>%
    left_join(isolates_traits) %>%
    filter(!genome_id %in% c("g20", "g28"))
list_traits <- c(paste0(rep(c("r", "lag", "maxOD"), each = 4), "_", rep(c(25, 30, 35, 40), 3), "c"),
                 "root_biomass_mg", "shoot_biomass_mg", "nodule_count")

# Core gene tree
list_scaled_branches <- c(37:39, 14:16, 50)
p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    geom_tiplab(align = T, hjust = -0.1, size = 3) +
    geom_tippoint(aes(color = site_group)) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(40, -5, 0) *1e-4, nudge_y = c(0, 1, -1), hjust = 1, fontface = "italic") +
    geom_treescale(x = 0, y = 28, width = 0.001) +
    scale_linetype_manual(values = c(1,5)) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(expand = c(0,0.001)) +
    theme_tree() +
    theme(
        legend.position = "top",
        legend.key.spacing.y = unit(1, "mm")
    ) +
    guides(linetype = "none", color = guide_legend(title = NULL, nrow = 2)) +
    labs(title = "core gene")

# trait heatmap
ordered_tips <- p1$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)
ii <- isolates %>%
    mutate(genome_id = factor(genome_id, ordered_tips)) %>%
    select(genome_id, starts_with("r"), starts_with("lag"), starts_with("maxOD"), shoot_biomass_mg, root_biomass_mg, nodule_count) %>%
    pivot_longer(-genome_id, names_to = "trait") %>%
    mutate(trait = factor(trait, list_traits)) %>%
    # Scale
    group_by(trait) %>%
    mutate(value = (value - min(value, na.rm = T)) / (max(value, na.rm = T) - min(value, na.rm = T))) %>%
    mutate(trait_type = case_when(
        str_detect(trait, "r_") ~ "r",
        str_detect(trait, "lag_") ~ "lag",
        str_detect(trait, "maxOD_") ~ "yield",
        T ~ "symbiosis"
    )) %>%
    mutate(trait_type = factor(trait_type, c("r", "lag", "yield", "symbiosis"))) %>%
    mutate(trait = str_remove(trait, "r_|lag_|maxOD_")) %>%
    mutate(trait = case_when(
        trait == "shoot_biomass_mg" ~ "W[s]",
        trait == "root_biomass_mg" ~ "W[r]",
        trait == "nodule_count" ~ "N",
        T ~ trait
    )) %>%
    mutate(trait = factor(trait, c("25c", "30c", "35c", "W[s]", "W[r]", "N"))) %>%
    ungroup()

p2 <- ii %>%
    filter(trait != "40c") %>%
    ggplot() +
    geom_tile(aes(x = trait, y = genome_id, fill = value)) +
    scale_fill_gradient2(low = "gold", mid = "#d55e00", high = "black", na.value = "snow", midpoint = 0.5, name = "scaled") +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    facet_nested(.~ trait_type, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.placement = "outside",
        panel.border = element_rect(color = "black", fill = NA),
        panel.spacing.x = unit(1, "mm"),
        legend.key.height = unit(3, units = "lines"),
        plot.margin = unit(c(10,0,0,0), "mm"),
    ) +
    labs()

# Phylogenetic signals
tmp <- traits_ps %>%
    filter(tree == "core") %>%
    mutate(trait = factor(trait, list_traits)) %>%
    mutate(trait_type = case_when(
        str_detect(trait, "r_") ~ "r",
        str_detect(trait, "lag_") ~ "lag",
        str_detect(trait, "maxOD_") ~ "yield",
        T ~ "symbiosis"
    )) %>%
    mutate(trait_type = factor(trait_type, c("r", "lag", "yield", "symbiosis"))) %>%
    mutate(trait = str_remove(trait, "r_|lag_|maxOD_")) %>%
    mutate(trait = case_when(
        trait == "shoot_biomass_mg" ~ "W[s]",
        trait == "root_biomass_mg" ~ "W[r]",
        trait == "nodule_count" ~ "N",
        T ~ trait
    )) %>%
    mutate(trait = factor(trait, c("25c", "30c", "35c", "W[s]", "W[r]", "N"))) %>%
    drop_na(trait) %>%
    ungroup()
temp_p <- tmp %>%
    select(trait_type, trait, k = p_k, lambda = p_lambda) %>%
    pivot_longer(cols = c(k, lambda)) %>%
    mutate(name = ifelse(name == "k", "K", "\u03bb")) %>%
    mutate(value = case_when(
        value < 0.001 ~ "***",
        value < 0.01 ~ "**",
        value < 0.05 ~ "*",
        value >= 0.05 ~ "n.s.",
    ))
p3 <- tmp %>%
    select(trait, trait_type, k, lambda) %>%
    pivot_longer(cols = c(k, lambda)) %>%
    mutate(name = ifelse(name == "k", "K", "\u03bb")) %>%
    ggplot() +
    geom_tile(aes(x = trait, y = name, fill = value)) +
    geom_text(data = temp_p, aes(x = trait, y = name, label = value), size = 3) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "snow", high = "#AE211F", na.value = "black", name = "scaled") +
    facet_nested(.~ trait_type, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 2),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.spacing.x = unit(1, "mm"),
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.key.height = unit(0.5, units = "lines")
    ) +
    guides() +
    labs()

# p <- plot_grid(p1, p2, axis = "tb", align = "h", labels = LETTERS[1:2], scale = 1, label_x = c(0, -0.05)) +
#     theme(plot.background = element_rect(color = NA, fill = "white"))

p4 <- ggplot() + theme_void() + get_plot_component(p1, "guide-box", return_all = TRUE)

pp <- align_plots(p1+guides(color = "none"), p2, align = 'h', axis = "tb", greedy = T)
p <- plot_grid(pp[[1]], pp[[2]], NULL, p3, axis = "l", align = "v", labels = c("A", "B", "", "C"), scale = 1,
               label_x = c(0, -0.05, 0, -0.05), rel_heights = c(1, 0.2), rel_widths = c(0.6, 1)) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
p <- p + draw_plot(p4, scale = 0, x = 0.2, y = 0.1, halign = 0, valign = 0)

ggsave(here::here("plots/Fig5.png"), p, width = 10, height = 6)




if (F) {
    nol <- tibble(node = c(32, 34, 35), node_sp = c("Ensifer meliloti", "Ensifer medicae", "Ensifer spp."))
    # Tree
    p1 <- tr_acce %>%
        as_tibble() %>%
        left_join(rename(isolates, label = genome_id)) %>%
        as.treedata() %>%
        ggtree() +
        geom_tiplab(hjust = 1, offset = 0.001, as_ylab = T) +
        geom_tippoint(aes(color = site_group)) +
        #geom_label_repel(aes(label = node)) +
        geom_label2(aes(subset=(node %in% c(32,34,35))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
        geom_label2(aes(subset=(node %in% c(32,34,35))), label = c("Ensifer meliloti", "Ensifer medicae", "Ensifer spp."), label.size = 0, fill = NA, nudge_x = -.02, nudge_y = 1, hjust = 1) +
        annotate("segment", x = 0.1, xend = 0.2, y = 3, yend = 3, arrow = arrow(angle = 90, ends = "both", length = unit(1, "mm"))) +
        annotate("text", x = 0.15, y = 2, label = "0.1") +
        # geom_label2(aes(subset=(node %in% c(32,34,35))), label = c("Ensifer meliloti", "Ensifer medicae", "Ensifer spp."),
        #             label.r = unit(0.3, "lines"), label.size = 0, fill = c("maroon", "steelblue", "grey"), alpha = 0.5) +
        #geom_hilight(data = nol, aes(node = node, fill = node_sp), inherit.aes = F) +
        #geom_label(data = nol, aes(x = node, fill = node_sp), inherit.aes = F) +
        #geom_point2(aes(subset=(node %in% c(32,34,35))), label = c(1,2,3), shape = 23, size = 5, fill='red') +
        #geom_nodepoint(aes(x = node)) +
        #geom_nodepoint(data = nol, aes(x = node, fill = node_sp), inherit.aes = F) +
        #geom_nodelab(aes(label = node)) +
        #scale_color_manual(values = species_colors) +
        #scale_fill_manual(values = species_colors) +
        scale_x_continuous(expand = c(0,.001)) +
        #scale_x_continuous(limits = c(0, 1)) +
        theme_tree() +
        #theme_classic() +
        theme(
            legend.position = c(0.3, 1),
            legend.background = element_rect(color = "black", fill = "white"),
            axis.line.y = element_blank(),
            axis.ticks.length = unit(0, "mm"),
            axis.text.y = element_text(size = 10)
        ) +
        guides(color = guide_legend(title = NULL))
}
