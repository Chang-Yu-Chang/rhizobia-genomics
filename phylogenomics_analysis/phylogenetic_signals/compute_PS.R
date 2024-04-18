#' This scripts implements the phylogenetic comparative analysis

renv::load()
library(tidyverse)
library(cowplot)
library(phytools)
library(phangorn) # For rooting the tree
library(ggtree)
library(tidytree)
library(geiger) # For checking names
#library(nlme) # For PGLS
source(here::here("metadata.R"))

# Traits
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_traits <- read_csv(paste0(folder_data, "phenotypes_analysis/isolates_traits.csv"))
isolates <- isolates %>%
    left_join(isolates_contigs) %>%
    left_join(isolates_traits) %>%
    filter(!genome_id %in% c("g20", "g28"))

list_traits <- c(paste0(rep(c("r", "lag", "maxOD"), each = 4), "_", rep(c(25, 30, 35, 40), 3), "c"),
                 "root_biomass_mg", "shoot_biomass_mg", "nodule_count")

# Tree
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
#tr, tr_acce, gpatl
# tr <- read.tree(paste0(folder_data, "genomics/mltree/isolates_core_b/aln.treefile"))
# list_others <- c(paste0("g", c(20, 28, 38:43)), "em1022", "usda1106", "em1021", "wsm419")
# tr <- tr %>% drop.tip(list_others)

## Compute phylogenetic signals
list_traits <- str_subset(names(isolates), "0c$|5c$|nodule|biomass|contig_length")
compute_ps <- function(tree, genome_id, trait_value) {
    set.seed(1)
    temp <- setNames(trait_value, genome_id)
    ps1 <- phylosig(tree, temp, test = T, nsim = 1000, method = "K")
    ps2 <- phylosig(tree, temp, test = T, nsim = 1000, method = "lambda")
    return(tibble(k = ps1$K, p_k = ps1$P, lambda = ps2$lambda, p_lambda = ps2$P))
}

tb1 <- tibble(trait = list_traits) %>%
    rowwise() %>%
    mutate(result = list(compute_ps(tr, isolates$genome_id, isolates[[trait]]))) %>%
    unnest(result)

tb2 <- tibble(trait = list_traits) %>%
    rowwise() %>%
    mutate(result = list(compute_ps(tr_acce, isolates$genome_id, isolates[[trait]]))) %>%
    unnest(result)

traits_ps <- bind_rows(mutate(tb1, tree = "core"), mutate(tb2, tree = "gcv"))
write_csv(traits_ps, paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))

# 1. Plot tree with traits ----
traits_ps <- read_csv(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))
traits_ps <- traits_ps %>%
    mutate(across(c(k, p_k, lambda, p_lambda), function(x) round(x, 3)))


p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(color = species), align = T) +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = "top"
    )

ordered_tips <- p1$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)

p2 <- isolates %>%
    mutate(genome_id = factor(genome_id, ordered_tips)) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = r_30c)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2)
    ) +
    guides() +
    labs()

tmp <- traits_ps %>%
    mutate(across(c(k, p_k, lambda, p_lambda), function(x) round(x, 3))) %>%
    filter(trait == "r_30c", tree == "core")

p_ps <- tmp %>%
    ggplot() +
    geom_text(x = 0.5, y = 0.5, aes(label = paste("K = ", k, ", P = ", p_k, "\n\u03bb = ", lambda, ", P = ", p_lambda)), hjust = 0) +
    theme_classic()


p <- plot_grid(p1, p2, axis = "tb", align = "h")
ggsave(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/01-trait_ps_core.png"), p, width = 10, height = 6)

# 2. plot the GCV tree ----
traits_ps <- read_csv(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))

p1 <- tr_acce %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(color = species), align = T) +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = "top"
    )

ordered_tips <- p1$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)

p2 <- isolates %>%
    mutate(genome_id = factor(genome_id, ordered_tips)) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = r_25c)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2)
    ) +
    guides() +
    labs()

tmp <- traits_ps %>%
    mutate(across(c(k, p_k, lambda, p_lambda), function(x) round(x, 3))) %>%
    filter(trait == "r_25c", tree == "gcv")

p_ps <- tmp %>%
    ggplot() +
    geom_text(x = 0.5, y = 0.5, aes(label = paste("K = ", k, ", P = ", p_k, "\n\u03bb = ", lambda, ", P = ", p_lambda)), hjust = 0) +
    theme_classic()

p <- plot_grid(p1, p2, axis = "tb", align = "h")
ggsave(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/02-trait_ps_gcv.png"), p, width = 10, height = 6)

# 3. Plot the the core and GCV trees with heatmap  ----
## GCV tree
p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(align = T) +
    geom_nodelab(aes(label = label)) +
    #geom_nodelab(aes(label = node)) +
    #geom_highlight(node = 38, fill = species_colors["medicae"]) +
    #geom_tippoint(aes(color = species)) +
    #scale_color_manual(values = species_colors) +
    scale_x_continuous(limits = c(0, 0.65)) +
    #theme_classic() +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides() +
    labs()


## traits heatmap
ordered_tips <- p1$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)

p2 <- isolates %>%
    mutate(genome_id = factor(genome_id, ordered_tips)) %>%
    select(genome_id, starts_with("r"), starts_with("lag"), starts_with("maxOD"), shoot_biomass_mg, root_biomass_mg, nodule_count) %>%
    pivot_longer(-genome_id, names_to = "trait") %>%
    mutate(trait = factor(trait, list_traits)) %>%
    # Scale
    group_by(trait) %>%
    mutate(value = (value - min(value, na.rm = T)) / (max(value, na.rm = T) - min(value, na.rm = T))) %>%
    ggplot() +
    geom_tile(aes(x = trait, y = genome_id, fill = value)) +
    scale_fill_gradient2(low = "steelblue", mid = "snow", high = "maroon", na.value = "black", midpoint = 0.5, name = "scaled") +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.spacing = unit(0, "mm"),
        legend.key.height = unit(3, units = "lines"),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    #guides(fill = guide_legend(title = "scaled value")) +
    labs()

## phylogenetic signals
tmp <- traits_ps %>%
    filter(tree == "gcv") %>%
    mutate(trait = factor(trait, list_traits)) %>%
    drop_na(trait)
temp_p <- tmp %>%
    select(trait, k = p_k, lambda = p_lambda) %>%
    pivot_longer(cols = -trait) %>%
    mutate(name = ifelse(name == "k", "K", "\u03bb")) %>%
    mutate(value = case_when(
        value < 0.001 ~ "***",
        value < 0.01 ~ "**",
        value < 0.05 ~ "*",
        value >= 0.05 ~ "n.s.",
    ))

p3 <- tmp %>%
    select(trait, k, lambda) %>%
    pivot_longer(cols = c(k, lambda)) %>%
    mutate(name = ifelse(name == "k", "K", "\u03bb")) %>%
    ggplot() +
    geom_tile(aes(x = trait, y = name, fill = value)) +
    geom_text(data = temp_p, aes(x = trait, y = name, label = value)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient2(low = "steelblue", mid = "snow", high = "maroon", na.value = "black", midpoint = 0.5, name = "scaled") +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.spacing = unit(0, "mm"),
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.key.height = unit(3, units = "lines")
    ) +
    guides(fill = "none") +
    labs()

p <- plot_grid(p1, p2, NULL, p3, nrow = 2, axis = "", align = "vh", rel_heights = c(1, 0.3), scale =0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/03-traits_ps_core.png"), p, width = 10, height = 8)



