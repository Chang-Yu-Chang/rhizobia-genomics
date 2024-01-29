#' This script plot the NJ trees

renv::load()
library(tidyverse)
library(janitor)
library(ape)
library(ggtree)
library(tidytree)
source(here::here("analysis/00-metadata.R"))

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
#tree_symbiosis <- dists %>% select(genome_id1, genome_id2, d_symbiosis) %>% drop_na() %>% make_tree()
#save(tree_ani, tree_kmer, tree_growth, file = paste0(folder_data, "temp/32-trees.RData"))

# 1. Plot the kmer tree 
p <- tree_kmer %>%
    left_join(rename(genomes, label = genome_id), by = "label") %>%
    as.treedata() %>%
    ggtree() + # %<+% contigs_label +
    geom_tiplab(aes(color = batch_name)) +
#    scale_color_manual(values = rhizobia_population_colors) +
#    scale_x_continuous(limits = c(-1.1, 0.3)) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    #guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "contig gene content")
p
ggsave(paste0(folder_data, "temp/32a-03-tree_kmer.png"), p, width = 8, height = 4)
