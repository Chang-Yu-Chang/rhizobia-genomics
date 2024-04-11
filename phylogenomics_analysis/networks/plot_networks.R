#' This script plots the network based one distances

renv::load()
library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(ggforce)
library(graphlayouts) # for some graph layout algorithms that are not available in igraph.
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))
dist_genetics <- read_csv(paste0(folder_data, "genomics_analysis/dist_genetics.csv"))

# Create networks
nodes <- isolates_contigs %>% select(genome_id, everything()) %>%
    filter(genome_id %in% isolates$genome_id) %>%
    mutate(species = factor(species, unique(species))) %>%
    arrange(species) %>%
    filter(!genome_id %in% c("g20", "g28"))
edges <- dist_genetics %>%
    filter(genome_id1 %in% isolates$genome_id, genome_id2 %in% isolates$genome_id) %>%
    rename(from = genome_id1, to = genome_id2) %>%
    filter(!from %in% c("g20", "g28"), !to %in% c("g20", "g28"))
g <- tbl_graph(nodes, edges, directed = F)

save(g, file = paste0(folder_data, "phylogenomics_analysis/networks/networks.rdata"))

# 1. Plot kmer networks
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
p <- edges %>%
    mutate(to = factor(to, nodes$genome_id)) %>%
    mutate(from = factor(from, rev(nodes$genome_id))) %>%
    ggplot() +
    geom_tile(aes(x = to, y = from, fill = d_kmer)) +
    scale_fill_viridis(alpha = 0.8, direction = -1, begin = 0, end = 1, breaks = seq(0,1,0.25)) +
    scale_x_discrete(position = "top", drop = F) +
    scale_y_discrete(drop = F) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text = element_blank()
    ) +
    labs(title = "Jaccard distance")

ggsave(paste0(folder_data, "phylogenomics_analysis/networks/01-kmers_heatmap.png"), p, width = 4, height = 3)

# 2. Plot kmer networks
hist(dist_genetics$d_kmer)
thr_kmer <- 0.85
p <- g %>%
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
ggsave(paste0(folder_data, "phylogenomics_analysis/networks/02-kmers_network.png"), p, width = 4, height = 3)
























