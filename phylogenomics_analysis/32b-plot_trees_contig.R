#' This script plots the NJ trees at the contig level

renv::load()
library(tidyverse)
library(cowplot)
library(ggsci)
library(ape)
library(ggtree)
library(tidytree)
source(here::here("analysis/00-metadata.R"))

load(paste0(folder_data, "temp/32-trees.RData"))
contigs <- read_csv(paste0(folder_data, "temp/14-contigs.csv"))

# 1. Plot the contig kmer tree
plot_tree1 <- function (tree, gtitle) {
    tree %>%
        as_tibble() %>%
        left_join(rename(contigs, label = contig_id), by = "label") %>%
        as.treedata() %>%
        ggtree(layout = "circular") +
        geom_tiplab(aes(color = replicon), size = 5) +
        theme_tree(legend.position = 'centre') +
        theme(
            legend.position = c(0.1,0.95),
            legend.background = element_rect(fill = NA, color = NA),
            plot.background = element_rect(fill = NA)
        ) +
        guides(color = guide_legend(title = NULL)) +
        labs(title = gtitle)
}

p1 <- list_trees_contigs$meliloti_kmer %>% plot_tree1("kmer meliloti")
p2 <- list_trees_contigs$medicae_kmer %>% plot_tree1("kmer medicae")

p <- plot_grid(p1, p2,
    nrow = 1, scale = 0.95) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_data, "temp/32b-01-trees_kmer.png"), p, width = 20, height = 10)

# 2. Plot the contig jaccard tree
p1 <- list_trees_contigs$meliloti_jaccard %>% plot_tree1("jaccard meliloti")
p2 <- list_trees_contigs$medicae_jaccard %>% plot_tree1("jaccard medicae")

p <- plot_grid(p1, p2,
    nrow = 1, scale = 0.95) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_data, "temp/32b-02-trees_jaccard.png"), p, width = 20, height = 10)

