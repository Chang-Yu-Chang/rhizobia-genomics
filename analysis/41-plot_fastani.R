#' Plot the genomic distance

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 0. read data
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    mutate(genome_id = paste0("g", ID))

tb_ani <- read_table(paste0(folder_data, "temp/plasmidsaurus/summary/34-medaka/ani.out"), show_col_types = F,
                     col_names = c("g_A", "g_B", "ani", "aligned_match", "total_sequence")) %>%
    mutate(g_A = str_remove(g_A, "/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/34-medaka/")) %>%
    mutate(g_B = str_remove(g_B, "/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/34-medaka/")) %>%
    mutate(g_A = str_remove(g_A, ".fasta")) %>%
    mutate(g_B = str_remove(g_B, ".fasta")) %>%
    mutate(g_A = str_replace(g_A, "consensus_g", "Chang_Q5C_")) %>%
    mutate(g_B = str_replace(g_B, "consensus_g", "Chang_Q5C_")) %>%
    mutate(g_A = factor(g_A, c(paste0("Chang_Q5C_", 1:20), "em1021", "em1022", "wsm419"))) %>%
    mutate(g_B = factor(g_B, c(paste0("Chang_Q5C_", 1:20), "em1021", "em1022", "wsm419")))


# Clean up
isolates_RDP %>%
    select(genome_id, Genus)

# Write the ani for only the Ensifers, including the three ncbi strains
tb_ani_ensifer <- tb_ani %>%
    mutate(ani = ani/100) %>%
    filter(!g_A %in% paste0("Chang_Q5C_", c(1, 7, 12, 14, 18))) %>%
    filter(!g_B %in% paste0("Chang_Q5C_", c(1, 7, 12, 14, 18))) %>%
    select(g_A, g_B, ani) %>%
    arrange(g_A, g_B) %>%
    pivot_wider(names_from = g_A, values_from = ani) %>%
    rename(key = g_B)
write_delim(tb_ani_ensifer, file = paste0(folder_data, "temp/plasmidsaurus/summary/34-medaka/ani_ensifer.txt"), delim = "\t")
#read_table(paste0(folder_data, "temp/plasmidsaurus/summary/34-medaka/ani_ensifer.txt"))

# 1. plot the pairwise ani estimate
p <- tb_ani %>%
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

ggsave(paste0(folder_data, "temp/41-01-ani_heatmap.png"), p, width = 6, height = 5)

# 2. plot the histogram of ani
p <- tb_ani %>%
    ggplot() +
    geom_histogram(aes(x = ani), color = "black", fill = "white", binwidth = 1) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/41-02-ani_histogram.png"), p, width = 4, height = 4)

# 2. plot the heatmap with labels


























