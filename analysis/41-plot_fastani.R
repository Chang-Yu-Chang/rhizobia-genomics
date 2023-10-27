#' Plot the genomic distance

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))


#
tb_abi <- read_table(paste0(folder_data, "temp/plasmidsaurus/summary/04-medaka/test.out"),
           col_names = c("g_A", "g_B", "ani", "aligned_match", "total_sequence")) %>%
    mutate(g_A = str_replace(g_A, "/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/04-medaka/consensus_", "")) %>%
    mutate(g_A = str_replace(g_A, ".fasta", "")) %>%
    mutate(g_B = str_replace(g_B, "/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/04-medaka/consensus_", "")) %>%
    mutate(g_B = str_replace(g_B, ".fasta", "")) %>%
    mutate(g_A = factor(g_A, paste0("g", 1:20))) %>%
    mutate(g_B = factor(g_B, paste0("g", 20:1)))


# 1. plot the pairwise ani estimate
p <- tb_abi %>%
    ggplot() +
    geom_tile(aes(x = g_A, y = g_B, fill = ani)) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "white", high = "maroon") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/41-01-pairwise_ani.png"), p, width = 6, height = 5)
