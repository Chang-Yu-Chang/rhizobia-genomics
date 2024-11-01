#' This script plots the consensus tree from iqtree

library(tidyverse)
library(cowplot)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))


# 1. Plot core gene tree ----
set_name = "elev_med"
tt <- read_gpas(set_name)

p1 <- tr_elev_med_core %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = population), hjust = 0, align = T) +
    #geom_nodelab(aes(label = label), color = "black", hjust = 0, size = 2) +
    geom_nodepoint(color = "grey70", shape = 16, alpha = .5, size = 3) +
    scale_color_manual(values = population_colors) +
    coord_cartesian(clip = "off") +
    geom_treescale() +
    theme_tree() +
    theme(
        legend.position = "none",
        legend.position.inside = c(0.1, 0.9),
        legend.title = element_blank(),
        plot.title = element_text(size = 8),
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides() +
    labs(title = paste0("Elevation S. medicae\n", nrow(tt$list_sccg), " single-copy core genes"))
p1
set_name = "urbn_mel"
tt <- read_gpas(set_name)

p2 <- tr_urbn_mel_core %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = population), hjust = 0, align = T) +
    #geom_nodelab(aes(label = label), color = "black", hjust = 0, size = 2) +
    geom_nodepoint(color = "grey70", shape = 16, alpha = .5, size = 3) +
    scale_color_manual(values = population_colors) +
    geom_treescale() +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        legend.position.inside = c(0.8, 0.3),
        legend.title = element_blank(),
        plot.title = element_text(size = 8),
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides() +
    labs(title = paste0("Urbanization S. meliloti\n", nrow(tt$list_sccg), " single-copy core genes"))

p <- plot_grid(p1, p2, nrow = 2)

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/01-combined_sccg.png"), p, width = 3, height = 6)





if (F) {

# 2. Plot Ensifer species tree ----
## meliloti tree
list_non_medicae <- isolates_contigs %>% filter(species != "medicae") %>% pull(genome_id)
list_non_meliloti <- isolates_contigs %>% filter(species != "meliloti") %>% pull(genome_id)

p1 <- tr_seq_core %>%
    drop.tip(list_non_meliloti) %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0, align = TRUE) +
    geom_nodelab(aes(label = label), color = "black", hjust = 0, size = 2) +
    geom_treescale(x = 0, y = 14, width = 0.001) +
    scale_color_manual(values = species_colors) +
    scale_x_continuous(expand = c(0,0.0005), limits = c(NA, .0075)) +
    theme_tree() +
    #theme_classic() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(color = "none") +
    labs(title = "E. meliloti")

## medicae tree
p2 <- tr_seq_core %>%
    drop.tip(list_non_medicae) %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0, align = TRUE) +
    geom_nodelab(aes(label = label), color = "black", hjust = 0, size = 2) +
    geom_treescale(x = 0, y = 11, width = 0.001) +
    scale_color_manual(values = species_colors) +
    scale_x_continuous(expand = c(0,0.0005), limits = c(NA, .0055)) +
    theme_tree() +
    #theme_classic() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(color = "none") +
    labs(title = "E. medicae")

p <- plot_grid(p1, p2, scale = 0.9) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/02-core_mltrees_mm.png"), p, width = 8, height = 4)
}

if (F) {

# 4. Plot gene content trees ----
p <- tr_gpa_genomes %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0) +
    scale_color_manual(values = species_colors) +
    #scale_x_continuous(limits = c(0, 1)) +
    theme_tree() +
    #theme_classic() +
    labs()
ggsave(paste0(folder_data, "phylogenomics_analysis/trees/04-pap_tree.png"), p, width = 6, height = 4)

# 5. Plot core and gene content trees ----
# Plot core tree
p1 <- tr_seq_core %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    labs()

# Plot accessory tree
p2 <- tr_gpa %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0) +
    scale_color_manual(values = species_colors)

d1 <- p1$data
d2 <- p2$data
## reverse x-axis and set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
dd <- bind_rows(d1, d2) %>% filter(isTip)

p <- p1 + geom_tree(data = d2) +
    ggnewscale::new_scale_fill() +
    geom_line(data = dd, aes(x, y, group = label), color = "grey90", linetype = 1, linewidth = .2) +
    geom_tippoint(aes(color = species), size = 2) +
    geom_tippoint(data = d2, aes(color = species), size = 2) +
    scale_x_continuous(limits = c(-0.5, 3)) +
    annotate("text", x = c(-0.3, 2.5), y = c(30,30), label = c("core genes", "gene content")) +
    theme_tree() +
    #theme_classic() +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/05-matched_tree.png"), p, width = 7, height = 3)

# 6. Plot the PAP heatmap, order by sites ----
ii <- isolates %>% select(genome_id, population, population, species)
p <- gpatl %>%
    left_join(ii) %>%
    mutate(genome_id = factor(genome_id, rev(ii$genome_id))) %>%
    ggplot() +
    geom_rect(data = tibble(population = names(population_colors[-5])), aes(fill = population), xmin = 0, xmax = 30000, ymin = 0, ymax = 100, alpha = 0.2) +
    geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = population_colors) +
    facet_grid(population~., scales = "free_y", space = "free_y") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        strip.background = element_rect(color = NA, fill = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome", title = "Gene presence/absence")

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/06-pap_heatmap_sites.png"), p, width = 7, height = 5)

# 7. Plot the core and gene content trees, color by sites ----
p1 <- tr_seq_core %>%
    as_tibble() %>%
    left_join(rename(ii, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tippoint(aes(color = population)) +
    scale_color_manual(values = population_colors, name = NULL) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    labs()

# Plot accessory tree
p2 <- tr_gpa %>%
    as_tibble() %>%
    left_join(rename(ii, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = population), hjust = 0) +
    scale_color_manual(values = population_colors)

d1 <- p1$data
d2 <- p2$data
## reverse x-axis and set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
dd <- bind_rows(d1, d2) %>% filter(isTip)

p <- p1 + geom_tree(data = d2) +
    ggnewscale::new_scale_fill() +
    geom_line(data = dd, aes(x, y, group = label), color = "grey90", linetype = 1, linewidth = .2) +
    geom_tippoint(aes(color = population), size = 2) +
    geom_tippoint(data = d2, aes(color = population), size = 2) +
    scale_x_continuous(limits = c(-0.5, 3)) +
    annotate("text", x = c(-0.3, 2.5), y = c(30,30), label = c("core genes", "gene content")) +
    theme_tree() +
    theme(
        legend.background = element_rect(color = "black", linewidth = .5)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/07-matched_tree_sites.png"), p, width = 7, height = 3)

# 8. Plot gene content trees ----
p <- tr_gpa %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0) +
    #geom_nodelab(aes(label = label)) +
    scale_color_manual(values = species_colors) +
    #scale_x_continuous(limits = c(0, 1)) +
    theme_tree() +
    #theme_classic() +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/08-gpa_mltree.png"), p, width = 7, height = 5)


# 9. Plot core gene tree with scaled branch length ----
list_scaled_branches <- c(37:39, 14:16, 50)
p <- tr_seq_core %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    #ggtree(branch.length = "none") +
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 16, color = 1, size = 3, alpha = 0.3) +
    geom_tiplab(align = T, hjust = -0.1, size = 3) +
    #geom_nodelab(aes(label = node)) +
    #geom_tippoint(aes(color = species)) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(20, -5, 0) *1e-4, nudge_y = c(1, 1, -1), hjust = 1, fontface = "italic") +
    geom_treescale(x = 0, y = 28, width = 0.001) +
    #scale_color_manual(values = species_colors) +
    scale_linetype_manual(values = c(1,5)) +
    #scale_color_manual(values=c("black", "firebrick", "steelblue", 1,2,3))
    scale_x_continuous(expand = c(0,0.001)) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(linetype = "none") +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/09-core_mltree_scaled.png"), p, width = 6, height = 4)

# 10. Plot scaled gpa tree with scaled branch length ----
list_scaled_branches <- c(1,2,31,33,34,35,46)
p <- tr_gpa %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    #ggtree() +
    #ggtree(branch.length = "none") +
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 16, color = 1, size = 3, alpha = 0.3) +
    geom_tiplab(align = T, size = 3) +
    #geom_tippoint(aes(color = species)) +
    #geom_label(aes(x = branch, label = node)) +
    #geom_nodelab(aes(label = node)) +
    geom_label2(aes(subset=(node %in% c(32,35,46))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(aes(subset=(node %in% c(32,35,46))), label = c("Ensifer spp.", "Ensifer medicae", "Ensifer meliloti"), label.size = 0, fill = NA, nudge_x = c(60, -10, -10)*1e-3, nudge_y = c(1, 1, 1), hjust = 1, fontface = "italic") +
    geom_treescale(x = 0, y = 28, width = 0.01) +
    scale_color_manual(values = species_colors) +
    scale_linetype_manual(values = c(1,5)) +
    #scale_color_manual(values=c("black", "firebrick", "steelblue", 1,2,3))
    #scale_x_continuous(expand = c(0,0.001)) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(linetype = "none") +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/10-gpa_mltree_scaled.png"), p, width = 6, height = 4)

# 11. plot matched scaled trees ----
list_scaled_branches <- c(37:39, 14:16, 50)
scaling_factor = 20
pt1 <- tr_seq_core %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    # scaling factor for matching tree x scales
    mutate(branch.length = branch.length * scaling_factor) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    #geom_tiplab(align = T, hjust = -0.1, size = 3) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(-5, -5, 0) *1e-4*scaling_factor, nudge_y = c(1, 1, -1), hjust = 1, fontface = "italic") +
    #geom_treescale(x = 0, y = 25, width = 0.01) +
    scale_linetype_manual(values = c(1,5)) +
    scale_x_continuous(expand = c(0,0.001)) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(linetype = "none") +
    labs()

# Plot accessory tree
list_scaled_branches <- c(1,2,31,33,34,35,46)
pt2 <- tr_gpa %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    #geom_tiplab(align = T, size = 3) +
    geom_label2(aes(subset=(node %in% c(32,35,46))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(aes(subset=(node %in% c(32,35,46))), label = c("Ensifer spp.", "Ensifer medicae", "Ensifer meliloti"), label.size = 0, fill = NA, nudge_x = c(60, -10, -10)*1e-3, nudge_y = c(1, 1, 1), hjust = 1, fontface = "italic") +
    geom_treescale(x = 0, y = 28, width = 0.01) +
    scale_color_manual(values = species_colors) +
    scale_linetype_manual(values = c(1,5)) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(linetype = "none", color = "none") +
    labs()

d1 <- pt1$data
d2 <- pt2$data
## reverse x-axis and set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.05
dd <- bind_rows(d1, d2) %>% filter(isTip)

scale_width = 0.01
p <- pt1 +
    geom_tree(data = d2, aes(linetype = scaled_branch)) +
    ggnewscale::new_scale_fill() +
    geom_line(data = dd, aes(x, y, group = label), color = "grey90", linetype = 1, linewidth = .2) +
    # left tree
    geom_tippoint(aes(color = species), size = 2) +
    annotate("segment", x = -0.05, xend = -0.05 + scale_width, y = 25, yend = 25) +
    annotate("text", x = -0.05+scale_width/2, y = 25, label = format(scale_width / scaling_factor, scientific = F), vjust = -1, hjust = 0.5) +
    # right tree
    geom_tippoint(data = d2, aes(color = species), size = 2) +
    geom_nodepoint(data = d2, aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    geom_label2(data = d2, aes(subset=(node %in% c(32,35,46))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(data = d2, aes(subset=(node %in% c(32,35,46))), label = c("Ensifer spp.", "Ensifer medicae", "Ensifer meliloti"), label.size = 0, fill = NA, nudge_x = c(-10, 90, 80)*1e-3, nudge_y = c(1, 1, 1), hjust = 1, fontface = "italic") +
    scale_color_manual(values = species_colors) +
    annotate("segment", x = 0.5, xend = 0.5+scale_width, y = 25, yend = 25) +
    annotate("text", x = 0.5+scale_width/2, y = 25, label = scale_width, vjust = -1, hjust = 0.5) +
    scale_x_continuous(limits = c(-0.1, 0.55), expand = c(0,0)) +
    # Label
    annotate("text", x = c(-0.05, 0.5), y = c(30,30), label = c("core genes", "gene content")) +
    theme_tree() +
    theme(
        legend.background = element_rect(color = "black", linewidth = .5)
    ) +
    guides(color = "none") +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/11-matched_mltree_scaled.png"), p, width = 10, height = 4)

# 12. plot the gene content tree by replicon ----

plot_gpa_tree <- function (tr_gpa, replicon_type) {
    tr_gpa %>%
        as_tibble() %>%
        left_join(rename(contigs, label = contig_id)) %>%
        as.treedata() %>%
        ggtree() +
        geom_tiplab(aes(label = label, color = species), hjust = 0) +
        coord_cartesian(clip = "off") +
        scale_color_manual(values = species_colors) +
        theme_tree() +
        theme(
            legend.position = "top",
            plot.margin = margin(10, 50, 10, 10, unit = "mm")
        ) +
        labs(title = replicon_type)

}

p1 <- tr_gpa_genomes %>%
    drop.tip(c("g2", "g3", "g15")) %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0) +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = margin(10, 50, 10, 10, unit = "mm")
    ) +
    labs(title = "genome")

p2 <- plot_gpa_tree(tr_gpa_chrom %>% drop.tip(c("g2_contig_1", "g3_contig_2", "g15_contig_5")) , "chromosome")
p3 <- plot_gpa_tree(tr_gpa_psyma, "pSymA like")
p4 <- plot_gpa_tree(tr_gpa_psymb, "pSymB like")
p <- plot_grid(p1, p2, p3, p4, nrow = 2)

ggsave(paste0(folder_data, "phylogenomics_analysis/trees/12-gpa_mltrees.png"), p, width = 15, height = 10)



























}
