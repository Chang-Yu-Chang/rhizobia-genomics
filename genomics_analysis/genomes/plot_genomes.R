#' This script plots the the genomes by contigs

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv")) %>%
    left_join(isolates) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))
qcs <- read_csv(paste0(folder_data, "genomics_analysis/genomes/qcs.csv"))

# 1. Number of contigs
p <- genomes %>%
    #filter(genome_id != "g28") %>%
    group_by(genome_id) %>%
    summarize(n_contigs = n()) %>%
    ggplot() +
    geom_histogram(aes(x = n_contigs), color = "black", fill = "white", binwidth = 1) +
    scale_x_continuous(breaks = seq(0, 12, 1)) +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "# of contigs", title = "contigs > 10kb")

ggsave(paste0(folder_data, "genomics_analysis/genomes/01-n_contigs.png"), p, width = 4, height = 4)


# 2. plot the genome size by contigs
p <- genomes %>%
    group_by(genome_id) %>%
    mutate(contig_ordered = factor(1:n())) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = contig_length/10^6, group = contig_ordered, fill = contig_ordered), position = position_stack(reverse = T), color = "black", width = 0.7) +
    #scale_y_continuous(expand = c(0,0), limits = c(0, 8.5), breaks = 1:8) +
    scale_fill_manual(values = c("grey", rep("white", 20)) %>% setNames(1:21), breaks = c(1,2), labels = c("largest contig", "other contig")) +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "genome", y = "genome size (Mbp)", title = "contigs < 10kb removed")

ggsave(paste0(folder_data, "genomics_analysis/genomes/02-genome_size.png"), p, width = 15, height = 6)
