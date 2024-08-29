#' This script plots the genomic distance matrices and mantel test

renv::load()
library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
# library(ggforce)
# library(graphlayouts) # for some graph layout algorithms that are not available in igraph.
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
dist_genomes <- read_csv(paste0(folder_data, "genomics_analysis/distances/dist_genomes.csv"))
tb_mantel <- read_csv(paste0(folder_data, "genomics_analysis/distances/tb_mantel.csv"))
dist_genomes <- dist_genomes %>%
    mutate(genome_id1 = ordered(genome_id1, genomes$genome_id), genome_id2 = ordered(genome_id2, genomes$genome_id)) %>%
    filter(genome_id1 >= genome_id2) %>%
    mutate(genome_id1 = ordered(genome_id1, rev(genomes$genome_id)), genome_id2 = ordered(genome_id2, genomes$genome_id))

# Plot distance matrixes
plot_dist <- function(d_i) {
    dist_genomes %>%
        ggplot() +
        geom_tile(aes(x = genome_id2, genome_id1, fill = {{d_i}}), width = .8, height = .8) +
        geom_tile(aes(x = genome_id2, genome_id1), color = "black", fill = NA) +
        scale_x_discrete(position = "bottom", expand = c(0,0)) +
        scale_y_discrete(position = "left", expand = c(0,0)) +
        scale_fill_gradient(low = "snow", high = "darkred", name = "distance") +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "inside",
            legend.position.inside = c(.8, .8),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.margin = unit(c(0,0,0,10), "mm")
        ) +
        guides() +
        labs(x = "genome", y = "genome")
}

plist <- list(
    plot_dist(d_kmer) + ggtitle("k-mer (31-mer)"),
    plot_dist(d_sccg) + ggtitle("988 single-copy core genes"),
    plot_dist(d_ani) + ggtitle("1-ANI"),
    plot_dist(d_jaccard) + ggtitle("Jaccard distance on gene presence/absence")
)


# Plot the network for mantel test
p_mantel <- tb_mantel %>%
    mutate(corr = round(corr,2)) %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(node_label = str_replace(name, "d_", "")) %>%
    ggraph(layout = "linear", circular = T) +
    geom_edge_arc2(aes(edge_width = corr)) +
    #geom_node_point(size = 10, shape = 21, fill = "gold", stroke = 0) +
    #geom_node_label(aes(label = node_label), label.size = 0, fill = NA) +
    scale_edge_width(range = c(.5,2), name = "Mantel's r") +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
        plot.margin = unit(c(10,10,10,10), "mm"),
        panel.background = element_rect(fill = alpha("grey", 0.3), color = NA),
        legend.position = "inside",
        legend.position.inside = c(1,.5),
        legend.key.height = unit(0, "mm")
    ) +
    guides()


#
p_dist <- plot_grid(plotlist = plist, ncol = 2, align = "hv", axis = "tblr", labels = LETTERS[1:4], label_x = .2, label_y = .9, scale = 0.75) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

p <- ggdraw(p_dist) +
    draw_plot(p_mantel, x = .4, y = .4, width = .25, height = .25)

ggsave(paste0(folder_data, "genomics_analysis/distances/01-distance_matrices.png"), p, width = 10, height = 10)


 if (F) {

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
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

save(nodes, edges, g, file = paste0(folder_data, "phylogenomics_analysis/networks/networks.rdata"))

# 1. Plot kmer ranges ----
p <- dist_genetics %>%
    filter(genome_id1 != genome_id2) %>%
    mutate(across(starts_with("genome_id"), ordered)) %>%
    filter(genome_id1 < genome_id2) %>%
    ggplot() +
    geom_histogram(aes(x = d_kmer), fill = "white", color = "black") +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/networks/01-kmers_hist.png"), p, width = 4, height = 3)

# 2. Plot kmer networks ----
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
p <- edges %>%
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

ggsave(paste0(folder_data, "phylogenomics_analysis/networks/02-kmers_heatmap.png"), p, width = 4, height = 3)

# 3. Plot kmer networks ----
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
ggsave(paste0(folder_data, "phylogenomics_analysis/networks/03-kmer_network.png"), p, width = 4, height = 3)


# 4. plot ani histogram ----
p <- dist_genetics %>%
    filter(genome_id1 != genome_id2) %>%
    mutate(across(starts_with("genome_id"), ordered)) %>%
    filter(genome_id1 < genome_id2) %>%
    ggplot() +
    geom_histogram(aes(x = d_ani), fill = "white", color = "black") +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "phylogenomics_analysis/networks/04-ani_hist.png"), p, width = 4, height = 3)

# 5. Plot ani heatmap ----
p <- edges %>%
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

ggsave(paste0(folder_data, "phylogenomics_analysis/networks/05-ani_heatmap.png"), p, width = 4, height = 3)


# 6. Plot ani networks ----
thr_ani <- 0.1
p <- g %>%
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
ggsave(paste0(folder_data, "phylogenomics_analysis/networks/06-ani_network.png"), p, width = 4, height = 3)





















}

