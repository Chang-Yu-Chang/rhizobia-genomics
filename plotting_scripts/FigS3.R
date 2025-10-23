#' This script plots the figures

library(tidyverse)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
ani <- read_csv(paste0(folder_genomics, "taxonomy/ani.csv"))
load(paste0(folder_genomics, "pangenome/trees/trees.rdata"))
tbb <- read_csv(paste0(folder_genomics, "pangenome/tree_distance/tbb.csv"))

# Panel A. core gene ----

tr <- tbtr$tr[[1]]
edges_to_scale <- which(tr$edge[,2] %in% 44)
tr$edge.length[edges_to_scale] <- tr$edge.length[edges_to_scale]*0.1

p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(ani, label = genome_id)) %>%
    mutate(` ` = "") %>%
    mutate(highlight = ifelse(node %in% 44, T, F)) %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = organism_name), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = organism_name), shape = -1, size = -1) +
    geom_nodepoint(aes(color = highlight, size = highlight)) +
    scale_color_manual(values = c(species_colors, `TRUE` = "red")) +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 0)) +
    scale_x_continuous(limits = c(0, 0.045)) +
    geom_treescale(x = .001, y = 20) +
    facet_grid2(~` `) +
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
        plot.margin = unit(c(0,-3,0,0), "mm")
    ) +
    guides(color = "none", size = "none") +
    labs()

p1_1 <- isolates %>%
    left_join(select(ani, genome_id, organism_name)) %>%
    filter(genome_id %in% get_taxa_name(p1)) %>%
    select(genome_id, region, organism_name) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    ggplot() +
    geom_tile(aes(x = region, y = genome_id, fill = organism_name), color = "black", linewidth = .5) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors, breaks = c("Sinorhizobium meliloti", "Sinorhizobium medicae")) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.text = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.background = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
    labs()


# Panel B. gene content ----
p2 <- tbtr$tr[[2]] %>%
    as_tibble() %>%
    left_join(rename(ani, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = organism_name), hjust = 1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = organism_name), shape = -1, size = -1) +
    scale_x_reverse() +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = -18, y = 20) +
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

p2_1 <- isolates %>%
    left_join(select(ani, genome_id, organism_name)) %>%
    filter(genome_id %in% get_taxa_name(p1)) %>%
    select(genome_id, region, organism_name) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p2)))) %>%
    ggplot() +
    geom_tile(aes(x = region, y = genome_id, fill = organism_name), color = "black", linewidth = .5) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.background = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
    labs()

p2_2 <- isolates %>%
    filter(genome_id %in% get_taxa_name(p1)) %>%
    left_join(select(ani, genome_id, organism_name)) %>%
    select(genome_id, site, organism_name) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p2)))) %>%
    ggplot() +
    geom_tile(aes(x = "Site", y = genome_id), color = "black", fill = NA, linewidth = .5) +
    geom_text(aes(x = "Site", y = genome_id, label = site), size = 3) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        panel.grid = element_blank()
    ) +
    guides() +
    labs()

# Panel C. tree distance ----
tb_obv <- tbb %>%
    select(-rfds) %>%
    distinct(tr_type1, tr_type2, .keep_all = T) %>%
    filter(tr_type1 == "gpa", tr_type2 == "core")

p3 <- tbb %>%
    filter(tr_type1 == "gpa", tr_type2 == "core") %>%
    ggplot() +
    geom_histogram(aes(x = rfds), binwidth = .05, color = "black", fill = NA, boundary = 0) +
    geom_vline(data = tb_obv, aes(xintercept = rfd), color = "red", linetype = 2) +
    scale_x_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
    scale_y_continuous(position = "right") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = NA, color = NA)
    ) +
    guides() +
    labs(x = "Robinson–Foulds distances")


# Connecting lines ----
tips1 <- get_taxa_name(p1)
tips2 <- get_taxa_name(p2)
p4 <- tibble(genome_id = sort(get_taxa_name(p1))) %>%
    mutate(
        x1 = 1,
        x2 = 2,
        y1 = match(genome_id, rev(get_taxa_name(p1))),
        y2 = match(genome_id, rev(get_taxa_name(p2)))
    ) %>%
    ggplot() +
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), linewidth = .2, color = "grey10") +
    scale_y_continuous(expand = c(0,0), limits = c(0.5, 35.5)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme() +
    guides() +
    labs()



# ----
p <- plot_grid(
    p1, p1_1 + guides(fill = "none"),
    p4,
    p2_1 + guides(fill = "none"), p2_2, p2,
    nrow = 1,
    scale = .95, rel_widths = c(1,.1,.15,.1,.08,1),
    align = "h", axis = "tb",
    labels = c("A", "", "", "B")
) +
    draw_text("Single-copy core gene", x = .2, y = .95, size = 10, hjust = 0) +
    draw_text("Gene content variation", x = .63, y = .95, size = 10, hjust = 0) +
    draw_label("C", x = .79, y = .95, fontface = "bold") +
    draw_plot(get_legend(p1_1), x = -.4, y = .25) +
    draw_plot(p3, width = .2, height = .35, x = .8, y = .6) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/FigS3.png"), p, width = 10, height = 6)
