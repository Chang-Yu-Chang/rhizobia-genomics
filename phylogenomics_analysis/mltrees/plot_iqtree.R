#' This script plots the consensus tree from iqtree

renv::load()
library(tidyverse)
library(ape)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
tr <- read.tree(paste0(folder_data, "genomics/mltree/core_b/aln.treefile"))

# 1. Full tree
tr %>%
    ggtree() +
    geom_tiplab(aes(label = label), hjust = 0) +
    #geom_nodelab(aes(label = label), color = "red", hjust = 0) +
    theme_tree()


# 2. medicae tree
list_others <- c(paste0("g", 38:43), "em1022", "usda1106", "em1021", "wsm419")
list_non_medicae <- isolates_contigs %>% filter(species != "medicae") %>% pull(genome_id)
list_non_meliloti <- isolates_contigs %>% filter(species != "meliloti") %>% pull(genome_id)

tr %>%
    drop.tip(list_non_medicae) %>%
    drop.tip(list_others) %>%
    ggtree() +
    geom_tiplab(aes(label = label), hjust = 0) +
    geom_nodelab(aes(label = label), color = "red", hjust = 0) +
    theme_tree()

# 3. meliloti tree
tr %>%
    drop.tip(list_non_meliloti) %>%
    drop.tip(list_others) %>%
    ggtree() +
    geom_tiplab(aes(label = label), hjust = 0) +
    geom_nodelab(aes(label = label), color = "red", hjust = 0) +
    theme_tree()

# 4.
