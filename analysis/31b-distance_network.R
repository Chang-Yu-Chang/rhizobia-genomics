#' This script plots the network based one distances

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
library(ggsci)
library(tidygraph)
library(ggraph)
library(ggforce)
library(graphlayouts) # for some graph layout algorithms that are not available in igraph.
source(here::here("analysis/00-metadata.R"))

isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))
dists <- read_csv(paste0(folder_data, "temp/31-dists.csv")); dists <- dists %>% filter(genome_id1 != genome_id2)
dists_long <- read_csv(paste0(folder_data, "temp/31-dists_long.csv"))


# 1. Plot the ani network
nodes <- select(isolates_traits, genome_id, everything()) %>%
    left_join(isolates_contigs)
pairs_d <- dists %>%
    filter(genome_id1 %in% isolates_traits$genome_id, genome_id2 %in% isolates_traits$genome_id) %>%
    rename(from = genome_id1, to = genome_id2)

g <- tbl_graph(nodes, pairs_d, directed = F)
bb <- layout_as_backbone(g, keep = 0.2)

g %>%
    activate(edges) %>%
    mutate(ani_bin = ifelse(d_ani < 0.05, "similar", "nope")) %>%
    ggraph(x = bb$xy[, 1], y = bb$xy[, 2]) +
    geom_mark_hull(
        aes(x, y, group = species, fill = species),
        concavity = 4,
        expand = unit(2, "mm"),
        alpha = 0.25
    ) +
    geom_edge_link(aes(color = ani_bin)) +
    geom_node_point(aes(fill = species), shape = 21, color = "black", size = 5) +
    scale_edge_color_manual(values = c(similar = "black", nope = "NA")) +
    theme_graph() +
    labs()


pairs_d %>%
    left_join(select(nodes, from = genome_id, species1 = species)) %>%
    left_join(select(nodes, to = genome_id, species2 = species)) %>%
    mutate(p_species = ifelse(species1 == species2, "yes", "no")) %>%
    ggplot() +
    geom_histogram(aes(x = d_ani, fill = p_species)) +
    theme_classic() +
    theme() +
    guides() +
    labs()


# 2. Plot the kmer networks























