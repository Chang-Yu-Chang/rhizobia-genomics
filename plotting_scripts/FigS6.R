#' This script plots the kmer networks

renv::load()
library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(ggforce)
source(here::here("metadata.R"))

dist_genetics <- read_csv(paste0(folder_data, "genomics_analysis/dist_genetics.csv"))
load(file = paste0(folder_data, "phylogenomics_analysis/networks/networks.rdata"))

# Panel A. Plot kmer histogram
p1 <- dist_genetics %>%
    filter(genome_id1 != genome_id2) %>%
    mutate(across(starts_with("genome_id"), ordered)) %>%
    filter(genome_id1 < genome_id2) %>%
    ggplot() +
    geom_histogram(aes(x = d_kmer), fill = "white", color = "black") +
    theme_classic() +
    theme() +
    guides() +
    labs()

# Panel B. Plot kmer heatmap ----
thr_kmer <- 0.85
sort_pairs <- function (edges, upper = T) {
    tt <- edges %>%
        mutate(pair = 1:n()) %>%
        pivot_longer(cols = c(from, to)) %>%
        group_by(pair) %>%
        mutate(value = factor(value, nodes$genome_id)) %>%
        arrange(value)

    if (upper) {
        tt <- mutate(tt, name = c("from", "to"))
    } else tt <- mutate(tt, name = c("to", "from"))
    return(pivot_wider(tt, names_from = name, values_from = value))

}
p2 <- edges %>%
    mutate(to = factor(to, nodes$genome_id)) %>%
    mutate(from = factor(from, rev(nodes$genome_id))) %>%
    ggplot() +
    geom_tile(aes(x = to, y = from, fill = d_kmer)) +
    #scale_fill_viridis(alpha = 0.8, direction = -1, begin = 0, end = 1, breaks = seq(0,1,0.25)) +
    scale_fill_gradient2(low = "gold", mid = "palegreen4", high = "steelblue", midpoint = 0.8) +
    scale_x_discrete(position = "top", drop = F) +
    scale_y_discrete(drop = F) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text = element_blank()
    ) +
    labs(title = "Jaccard distance")

# Panel C. Plot kmer networks ----
p3 <- g %>%
    activate(edges) %>%
    mutate(kmer_bin = ifelse(d_kmer < thr_kmer, paste0("d_kmer<", thr_kmer), "nope")) %>%
    ggraph(layout = "linear", circular = T) +
    #geom_mark_hull(aes(x, y, group = species, fill = species), concavity = 4, expand = unit(3, "mm"), alpha = 0.25) +
    geom_edge_arc(aes(colour = kmer_bin), alpha = 0.1) +
    geom_node_point(aes(fill = species), shape = 21, color = "black", size = 5) +
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

# Panel D. plot ani histogram ----
p4 <- dist_genetics %>%
    filter(genome_id1 != genome_id2) %>%
    mutate(across(starts_with("genome_id"), ordered)) %>%
    filter(genome_id1 < genome_id2) %>%
    ggplot() +
    geom_histogram(aes(x = d_ani), fill = "white", color = "black") +
    theme_classic() +
    theme() +
    guides() +
    labs()

# Panel E. Plot ani heatmap ----
p5 <- edges %>%
    mutate(to = factor(to, nodes$genome_id)) %>%
    mutate(from = factor(from, rev(nodes$genome_id))) %>%
    ggplot() +
    geom_tile(aes(x = to, y = from, fill = d_ani)) +
    scale_fill_gradient2(low = "gold", mid = "palegreen4", high = "steelblue", midpoint = 0.1) +
    #scale_fill_viridis(alpha = 0.8, direction = -1, begin = 0, end = 0.3, breaks = seq(0,1,0.25)) +
    scale_x_discrete(position = "top", drop = F) +
    scale_y_discrete(drop = F) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text = element_blank()
    ) +
    labs(title = "ANI")

# Panel F. Plot ani networks ----
thr_ani <- 0.1
p6 <- g %>%
    activate(edges) %>%
    mutate(ani_bin = ifelse(d_ani < thr_ani, paste0("d_ani<", thr_ani), "nope")) %>%
    ggraph(layout = "linear", circular = T) +
    #geom_mark_hull(aes(x, y, group = species, fill = species), concavity = 4, expand = unit(3, "mm"), alpha = 0.25) +
    geom_edge_arc(aes(colour = ani_bin), alpha = 0.1) +
    geom_node_point(aes(fill = species), shape = 21, color = "black", size = 5) +
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

p <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, labels = LETTERS[c(1,3,5,2,4,6)], scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/FigS6.png"), p, width = 10, height = 5)


