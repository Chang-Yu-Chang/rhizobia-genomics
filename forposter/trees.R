#' This script plots the consensus tree from iqtree

renv::load()
library(tidyverse)
library(cowplot)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))


list_scaled_branches <- c(15,16,17,18,27,13,12,11,14)
tr_seq_core$tip.label
p1 <- tr_seq_core %>%
    drop.tip(isolates$genome_id[isolates$population == "PA"]) %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    left_join(rename(isolates, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, "scaled to 1%", "not scaled")) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch), layout="ellipse", branch.length="none", linewidth = 2) +
    #geom_tiplab(offset = .3, size = 8) +
    geom_tippoint(aes(color = site_group), size = 8) +
    geom_nodepoint(aes(label = highlight_boot, color = "bootstrap>95%"), shape = 16, size = 10, alpha = 0.3) +
    #geom_nodelab(aes(label = node)) +
    # geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    # geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(20, -5, 0) *1e-4, nudge_y = c(1, 1, -1), hjust = 1, fontface = "italic") +
    # geom_treescale(x = 0, y = 28, width = 0.001) +
    scale_color_manual(values = c(site_group_colors, `bootstrap>95%`="grey40"), name = NULL, breaks = c("high elevation", "low elevation", "bootstrap>95%")) +
    scale_linetype_manual(values = c(1,5), name = NULL) +
    #coord_cartesian(clip = "off") +
    coord_flip(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.2),
        legend.background = element_rect(color = NA, fill = NA),
        legend.box.background = element_rect(color = NA, fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        legend.spacing.y = unit(-3,"mm"),
        legend.text = element_text(size = 20),
        plot.margin = unit(c(0,10,0,0), "mm"),
        plot.background = element_blank()
    ) +
    guides(linetype = "none") +
    labs()

list_scaled_branches <- c(18, 31)
p2 <- tr_seq_core %>%
    drop.tip(isolates$genome_id[isolates$population == "VA"]) %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    left_join(rename(isolates, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, "scaled to 1%", "not scaled")) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch), layout="ellipse", branch.length="none", linewidth = 2) +
    #geom_nodepoint(aes(label = highlight_boot), shape = 16, color = 1, size = 3, alpha = 0.3) +
    #geom_tiplab(offset = .3, size = 8) +
    geom_tippoint(aes(color = site_group), size = 8) +
    geom_nodepoint(aes(label = highlight_boot, color = "bootstrap>95%"), shape = 16, size = 10, alpha = 0.3) +
    #geom_nodelab(aes(label = node)) +
    # geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    # geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(20, -5, 0) *1e-4, nudge_y = c(1, 1, -1), hjust = 1, fontface = "italic") +
    # geom_treescale(x = 0, y = 28, width = 0.001) +
    scale_color_manual(values = c(site_group_colors, `bootstrap>95%`="grey40"), name = NULL, breaks = c("suburban", "urban", "bootstrap>95%")) +
    scale_linetype_manual(values = c(1,5), name = NULL) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.9),
        legend.background = element_rect(color = NA, fill = NA),
        legend.box.background = element_rect(color = NA, fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        legend.spacing.y = unit(-3,"mm"),
        legend.text = element_text(size = 20),
        plot.margin = unit(c(0,10,0,0), "mm"),
        plot.background = element_blank()
    ) +
    guides(linetype = "none") +
    labs()
p <- plot_grid(p1, p2, nrow = 2, scale = .9)
ggsave(here::here("forposter/trees.pdf"), p, width = 8, height = 12)

