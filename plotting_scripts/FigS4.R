#'

library(tidyverse)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))

make_heatmap <- function (p1) {
    isolates %>%
        left_join(select(iso, genome_id, contig_species)) %>%
        select(genome_id, population, contig_species) %>%
        mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
        ggplot() +
        geom_tile(aes(x = population, y = genome_id, fill = contig_species), color = "black", linewidth = .5) +
        scale_x_discrete(expand = c(0,0), position = "top") +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_manual(values = species_colors, breaks = c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")) +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            legend.key.size = unit(3, "mm"),
            legend.key.spacing.y = unit(1, "mm"),
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            strip.placement = "outside",
            strip.clip = "off",
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
            panel.background = element_rect(color = "black", fill = NA),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(0,2,0,1), "mm")
        ) +
        guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
        labs()
}

# Panel A. core gene ----
nodes_to_scale <- c(38, 40, 1, 2, 41, 42, 54)
tr <- tbtr$tr[[1]]
tr <- root(tr, outgroup = "g2")
edges_to_scale <- which(tr$edge[,2] %in% nodes_to_scale)
tr$edge.length[edges_to_scale] <- tr$edge.length[edges_to_scale]*0.01

p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    mutate(highlight = ifelse(node %in% nodes_to_scale, T, F)) %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_color_manual(values = species_colors) +
    scale_x_continuous(limits = c(0, 0.0075)) +
    geom_treescale(x = .001, y = 15) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"),
        legend.position = "inside",
        legend.position.inside = c(.2, .8),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(color = "none") +
    labs()


# Panel B. gene content ----
p2 <- tbtr$tr[[2]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = 0, align = T, offset = 1, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 12, y = 30) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()

# Panel C. Structural variation ----
p3 <- tbtr$tr[[3]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = 0, align = T, offset = 1, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 12, y = 30) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()


# Panel D. ani ----
p4 <- tbtr$tr[[4]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = 0, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 0.02, y = 30) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()

# Panel E. kmer ----
p5 <- tbtr$tr[[5]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = 0, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 0.1, y = 30) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()


# ----
p <- plot_grid(
    p1 + ggtitle("single-copy core gene"), make_heatmap(p1) + guides(fill = "none"),
    p2 + ggtitle("Gene content variation"), make_heatmap(p2) + guides(fill = "none"),
    p3 + ggtitle("Structural variation"), make_heatmap(p3) + guides(fill = "none"),
    p4 + ggtitle("ANI"), make_heatmap(p4) + guides(fill = "none"),
    p5 + ggtitle("kmer"), make_heatmap(p5) + guides(fill = "none"),
    ncol = 6,
    scale = .9, rel_widths = c(1,.1,1,.1,1,.1),
    align = "h", axis = "tb",
    labels = c("A", "", "B", "", "C", "", "D", "", "E", "")
) +
    # draw_text("Single-copy core gene", x = .27, y = .95, size = 10, hjust = 0) +
    # draw_text("Gene content variation", x = .6, y = .95, size = 10, hjust = 0) +
    # draw_label("C", x = .79, y = .95, fontface = "bold") +
    # draw_plot(get_legend(p1_1), x = -.4, y = .25) +
    # draw_plot(p3, width = .2, height = .35, x = .8, y = .6) +
    theme(plot.background = element_rect(color = NA, fill = "white"))


ggsave(here::here("plots/FigS4.png"), p, width = 15, height = 10)

