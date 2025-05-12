#' This script plots the figures

library(tidyverse)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))

nodes_to_scale <- c(38, 40, 1, 2, 41, 42, 54)
tr <- tbtr$tr[[1]]
tr <- root(tr, outgroup = "g2")
edges_to_scale <- which(tr$edge[,2] %in% nodes_to_scale)
tr$edge.length[edges_to_scale] <- tr$edge.length[edges_to_scale]*0.01

tr %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    mutate(hilight = ifelse(node %in% nodes_to_scale, T, F)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = contig_species, fill = contig_species), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species, fill = contig_species), shape = -1, size = -1) +
    geom_nodepoint(aes(color = hilight), alpha = .5, size = 3) +
    #geom_nodelab(aes(label = node)) +
    scale_color_manual(values = species_colors) +
    # scale_color_manual(values = population_colors) +
    # scale_fill_manual(values = population_colors) +
    #scale_y_continuous(limits = c(1, length(tr_core$tip.label)), expand = c(0,.5)) +
    #coord_cartesian(clip = "off") +
    facet_grid2(~` `) +
    theme_tree() +
    theme(
        legend.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        plot.margin = unit(c(0,5,0,0), "mm")
    ) +
    #guides(fill = guide_legend(override.aes = list(label = "", color = NA, size = 2)), color = guide_legend(override.aes = list(size = 2, shape = 21))) +
    labs()


p <- plot_grid(
    p_tree1 + theme(legend.position = "none"),
    p_gpa1 + theme(legend.position = "none"),
    p_heat1,
    leg1, NULL, p_ngen1,
    NULL, NULL, NULL,
    p_tree2 + theme(legend.position = "none"),
    p_gpa2 + theme(legend.position = "none"),
    p_heat2,
    leg2, NULL, p_ngen2,
    ncol = 3,  align = "hv", axis = "blr", scale = .95,
    labels = c("A", "C", "E", rep("", 6), "B", "D", "F", rep("", 3)),
    rel_widths = c(1, 1, 2), rel_heights = c(10, 2, .5, 19, 2)
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig5.png"), p, width = 10, height = 6)
