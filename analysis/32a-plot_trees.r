#' This script plot the NJ trees

renv::load()
library(tidyverse)
library(cowplot)
library(ape)
library(ggtree)
library(tidytree)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
isolates_traits <- read_csv(paste0(folder_data, "temp/24-isolates_traits.csv"))

dist_genetics <- read_csv(paste0(folder_data, "temp/24-genetics.csv"))
dist_traits <- read_csv(paste0(folder_data, "temp/24-dist_traits.csv"))
dists <- read_csv(paste0(folder_data, 'temp/31-dists.csv'))
dists_long <- read_csv(paste0(folder_data, "temp/31-dists_long.csv"))


# Function for making tree
make_tree <- function(di) {
    #' This functions creates a phylo object from a long-format distance matrix
    colnames(di)[3] <- "dd"
    dist1 <- tibble(genome_id1 = di$genome_id1, genome_id2 = di$genome_id2, dd = di$dd)
    dist1_swapped <- tibble(genome_id1 = di$genome_id2, genome_id2 = di$genome_id1, dd = di$dd)
    dist1_swapped <- filter(dist1_swapped, genome_id1 != genome_id2)
    dist2 <- bind_rows(dist1, dist1_swapped)
    tree <- dist2 %>%
        pivot_wider(names_from = genome_id2, values_from = dd) %>%
        select(-genome_id1) %>%
        as.matrix() %>%
        nj() %>%
        as_tibble()
    return(tree)
}

tree_ani <- dists %>% select(genome_id1, genome_id2, d_ani) %>% make_tree()
tree_kmer <- dists %>% select(genome_id1, genome_id2, d_kmer) %>% make_tree()
tree_growth <- dists %>% select(genome_id1, genome_id2, d_growth) %>% drop_na() %>% make_tree()
tree_geo <- dists %>% select(genome_id1, genome_id2, d_geo) %>% drop_na() %>% make_tree()
#tree_symbiosis <- dists %>% select(genome_id1, genome_id2, d_symbiosis) %>% drop_na() %>% make_tree()
#save(tree_ani, tree_kmer, tree_growth, file = paste0(folder_data, "temp/32-trees.RData"))



# 1. Plot the trees by populations
plot_tree <- function (tree, gtitle) {
    tree %>%
        left_join(rename(isolates_traits, label = genome_id), by = "label") %>%
        as.treedata() %>%
        ggtree() +
        geom_tiplab(aes(color = population)) +
        geom_tippoint(aes(color = population)) +
        theme_tree(legend.position = 'centre') +
        theme(
            legend.position = c(0.2,0.9)
        ) +
        guides(color = guide_legend(override.aes = aes(label = ""), title = NULL)) +
        labs(title = gtitle)
}
p1 <- plot_tree(tree_ani, "ani")
p2 <- plot_tree(tree_kmer, "kmer")
p3 <- plot_tree(tree_growth, "growth")
p4 <- plot_tree(tree_geo, "geo")
p <- plot_grid(p1, p2, p3, p4, nrow = 2, scale = 0.9) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32a-01-trees.png"), p, width = 10, height = 10)

