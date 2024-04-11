#' This script plots the network based one distances

renv::load()
library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(ggforce)
library(ape)
library(tidytree)
library(ggtree)
library(proxy) # For computing jaccard distance
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
load(file = paste0(folder_data, "phylogenomics_analysis/networks/networks.rdata"))
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

# Panel A. kmer networks
thr_kmer <- 0.85
p1 <- g %>%
    activate(edges) %>%
    mutate(kmer_bin = ifelse(d_kmer < thr_kmer, paste0("d_kmer<", thr_kmer), "nope")) %>%
    ggraph(layout = "linear", circular = T) +
    #geom_mark_hull(aes(x, y, group = species, fill = species), concavity = 4, expand = unit(3, "mm"), alpha = 0.25) +
    geom_edge_arc(aes(colour = kmer_bin), alpha = 0.1) +
    geom_node_point(aes(fill = species), shape = 21, color = "black", size = 3) +
    scale_edge_color_manual(values = c("black", "NA")) +
    scale_fill_manual(values = species_colors) +
    theme_void() +
    theme(
        plot.background = element_rect(color = NA, fill = "white"),
        plot.margin = unit(c(1,1,1,1), "mm"),
        legend.position = "right"
    ) +
    guides(edge_colour = "none") +
    labs()

# Panel B. PAP heatmap
p2 <- gpatl %>%
    ggplot() +
    geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 5, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome")


# Panel C. matched tree
# Plot core tree
p3_1 <- tr %>%
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
p3_2 <- tr_acce %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0) +
    scale_color_manual(values = species_colors)

d1 <- p3_1$data
d2 <- p3_2$data
## reverse x-axis and set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
dd <- bind_rows(d1, d2) %>% filter(isTip)

p3 <- p3_1 + geom_tree(data = d2) +
    ggnewscale::new_scale_fill() +
    geom_line(data = dd, aes(x, y, group = label), color = "grey90", linetype = 1, linewidth = .2) +
    geom_tippoint(aes(color = species), size = 2) +
    geom_tippoint(data = d2, aes(color = species), size = 2) +
    scale_x_continuous(limits = c(-0.5, 3)) +
    annotate("text", x = c(-0.3, 2.5), y = c(30,30), label = c("core genes", "gene content")) +
    theme_tree() +
    #theme_classic() +
    labs()

p_top <- plot_grid(p1 + ggtitle("Jaccard distance"), p2 + ggtitle("Gene presence/absence"), nrow = 1, labels = c("A", "B"), scale = 0.8)
p <- plot_grid(p_top, p3, nrow = 2, labels = c("", "C"), scale = c(1, 0.9), rel_heights = c(1, 1.2)) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig4.png"), p, width = 8, height = 6)




