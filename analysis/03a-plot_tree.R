#' This script plots the tree based on the 16S multi sequence alignemnt

library(tidyverse)
library(ggtree)
library(tidytree)
library(tools)
source(here::here("analysis/00-metadata.R"))
options(ignore.negative.edge=TRUE) #The tree contained negative edge lengths. If you want to ignore the edges, you can set `options(ignore.negative.edge=TRUE)`, then re-run ggtree.

isolates_meta <- read_csv(paste0(folder_data, "temp/03-isolates_meta.csv"), show_col_types = F)
load(file = paste0(folder_data, "temp/03-tree_NJ.Rdata"))

#
isolates_meta <- isolates_meta %>%
    mutate(label = as.character(label))

# 1. Simple tree, label colored by location  ----
p <- tree_NJ %>%
    as_tibble() %>%
    left_join(isolates_meta) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(hjust = -0.2) +
    geom_tippoint(aes(color = Site), size = 2)+
    scale_color_manual(values = rhizobia_site_colors) +
    theme() +
    labs()
ggsave(paste0(folder_data, "temp/03a-01-tree_NJ.png"), p, width = 6, height = 4)

# 2. tree labeled colored by site, show genus
p <- tree_NJ %>%
    as_tibble() %>%
    left_join(isolates_meta) %>%
    as.treedata() %>%
    ggtree() +
    geom_label(aes(label = Genus, fill = Site), alpha = 0.5, hjust = 0) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_x_continuous(expand = c(0.2, 0)) +
    theme() +
    guides(label = "none") +
    labs()
ggsave(paste0(folder_data, "temp/03a-02-tree_NJ_genus.png"), p, width = 6, height = 5)
