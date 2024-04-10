#' This script plots the consensus tree from iqtree

renv::load()
library(tidyverse)
library(cowplot)
library(ape)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
tr <- read.tree(paste0(folder_data, "genomics/mltree/core_b/aln.treefile"))

# 1. Plot trees
## Full tree
p1 <- tr %>%
    drop.tip(list_others) %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    #ggtree(branch.length layout = "circular") +
    ggtree(layout = "roundrect") +
    geom_tiplab(aes(label = label, color = species), hjust = 0, align = TRUE) +
    geom_nodelab(aes(label = label), color = "maroon", hjust = 0) +
    #geom_nodelab(aes(label = node), color = "grey", hjust = 0) +
    # geom_hilight(node = 54, fill = "steelblue") +
    # geom_hilight(node = 52, fill = "maroon") +
    # geom_hilight(node = 37, fill = "pink") +
    #geom_hilight(aes(node = ), inherit.aes = FALSE) +
    #scale_color_continuous(low='white', high='red') +
    theme_tree()



## meliloti tree
list_others <- c(paste0("g", 38:43), "em1022", "usda1106", "em1021", "wsm419")
list_non_medicae <- isolates_contigs %>% filter(species != "medicae") %>% pull(genome_id)
list_non_meliloti <- isolates_contigs %>% filter(species != "meliloti") %>% pull(genome_id)

p2 <- tr %>%
    drop.tip(list_non_meliloti) %>%
    drop.tip(list_others) %>%
    #ggtree(layout = "slanted") +
    ggtree() +
    geom_tiplab(aes(label = label), hjust = 0, align = TRUE) +
    geom_nodelab(aes(label = label), color = "maroon", hjust = 0.2) +
    theme_tree()

## medicae tree
p3 <- tr %>%
    drop.tip(list_non_medicae) %>%
    drop.tip(list_others) %>%
    #ggtree(layout = "slanted") +
    ggtree() +
    geom_tiplab(aes(label = label), hjust = 0, align = TRUE) +
    geom_nodelab(aes(label = label), color = "maroon", hjust = 0.2) +
    theme_tree()

p <- plot_grid(p1, p2, p3, scale = 0.9)

ggsave(paste0(folder_data, "phylogenomics_analysis/mltrees/01-trees.png"), p, width = 10, height = 10)


# 4.













