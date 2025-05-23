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
tb2 <- read_csv(paste0(folder_data, "phylogenomics_analysis/tree_distance/tb2.csv"))

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
    guides(color = "none") +
    labs()

p1_1 <- isolates %>%
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
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
    labs()


# Panel B. gene content ----
p2 <- tbtr$tr[[2]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = 1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_x_reverse() +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = -18, y = 22) +
    coord_cartesian(clip = "off") +
    #theme_bw() +
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
    left_join(select(iso, genome_id, contig_species)) %>%
    select(genome_id, population, contig_species) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p2)))) %>%
    ggplot() +
    geom_tile(aes(x = population, y = genome_id, fill = contig_species), color = "black", linewidth = .5) +
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
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
    labs()

# Panel C. tree distance ----
tb_obv <- tb2 %>%
    select(-rfds) %>%
    distinct(tr_type1, tr_type2, .keep_all = T) %>%
    filter(tr_type1 == "gpa", tr_type2 == "core")

p3 <- tb2 %>%
    filter(tr_type1 == "gpa", tr_type2 == "core") %>%
    ggplot() +
    geom_histogram(aes(x = rfds), binwidth = .05, color = "black", fill = NA, boundary = 0) +
    geom_vline(data = tb_obv, aes(xintercept = rfd), color = "red", linetype = 2) +
    #facet_grid2(tr_type1 ~ tr_type2, switch = "y", render_empty = F) +
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

p4 <- tibble(genome_id = iso$genome_id) %>%
    mutate(
        x1 = 1,
        x2 = 2,
        y1 = match(genome_id, rev(get_taxa_name(p1))),
        y2 = match(genome_id, rev(get_taxa_name(p2)))
    ) %>%
    ggplot() +
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), linewidth = .2, color = "grey10") +
    scale_y_continuous(expand = c(0,0), limits = c(0.5, 38.5)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme() +
    guides() +
    labs()



# ----
p <- plot_grid(
    p1, p1_1 + guides(fill = "none"),
    p4,
    p2_1 + guides(fill = "none"), p2,
    nrow = 1,
    scale = .95, rel_widths = c(1,.1,.15,.1,1),
    align = "h", axis = "tb",
    labels = c("A", "", "B", "")
) +
    draw_text("Single-copy core gene", x = .27, y = .95, size = 10, hjust = 0) +
    draw_text("Gene content variation", x = .6, y = .95, size = 10, hjust = 0) +
    draw_label("C", x = .8, y = .97, fontface = "bold") +
    draw_plot(get_legend(p1_1), x = -.4, y = .25) +
    draw_plot(p3, width = .2, height = .35, x = .8, y = .6) +
    theme(plot.background = element_rect(color = NA, fill = "white"))


ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 6)


#
tt <- read_gpas()
nrow(tt$gpa) # 26544

core <- tt$gpatl %>%
    group_by(gene) %>%
    filter(value == 1) %>%
    count() %>%
    ungroup() %>%
    filter(n == max(n))
nrow(core) / nrow(tt$gpa) *100

tt$gpacl %>%
    filter(str_detect(gene, "nod")) %>%
    filter(genome_id %in% paste0("g", c(2,3,15)))

