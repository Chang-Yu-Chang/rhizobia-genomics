#' This script is to have a overview on the raw read statistics

library(tidyverse)
library(cowplot)
library(seqinr)
source(here::here("analysis/00-metadata.R"))

# 0. read the raw read data
#list_g <- paste0("Chang_Q5C_", c(2:6, 8:11, 13, 15:17, 19))
list_g <- paste0("Chang_Q5C_", 1:19)

list_g_contigs <- rep(list(NA), length(list_g))
for (i in 1:length(list_g)) {
    i=1
    fa <- read.fasta(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/04-medaka/consensus.fasta"))
    fa_len <- sapply(fa, length)
    list_g_contigs[[i]] <- tibble(genome_id = paste0("g", i), contig = names(fa_len), contig_length = fa_len)
}

g_contigs <- bind_rows(list_g_contigs) %>%
    arrange(genome_id, desc(contig_length)) %>%
    mutate(genome_id = factor(genome_id, paste0("g", 1:19)))

write_csv(g_contigs, paste0(folder_data, "temp/34-g_contigs.csv"))

# Some numbers
g_contigs %>%
    group_by(genome_id) %>%
    summarize(count = n())

# 1. plot the genome size by
p <- g_contigs %>%
    group_by(genome_id) %>%
    mutate(contig_ordered = factor(1:n())) %>%
    #summarize(genome_size = sum(contig_length)) %>%
    ggplot() +
    #geom_col(aes(x = genome_id, y = genome_size/10^6), color = "black", fill = "white") +
    geom_col(aes(x = genome_id, y = contig_length/10^6, group = contig_ordered, fill = contig_ordered), position = position_stack(reverse = T), color = "black") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 8.5), breaks = 1:8) +
    scale_fill_manual(values = c("grey", rep("white", 20)) %>% setNames(1:21), breaks = c(1,2), labels = c("largest contig", "other contig")) +
    theme_bw() +
    theme() +
    guides(fill = guide_legend(title = "")) +
    labs(x = "genome", y = "genome size (Mbp)")

ggsave(paste0(folder_data, "temp/34-01-genome_size.png"), plot = p, width = 8, height = 5)

# 2. plot the contig size for all genomes
p <- g_contigs %>%
    ggplot() +
    geom_histogram(aes(x = contig_length/10^6), color = "black", fill = "white") +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "contig size (Mbp)", y = "count")

ggsave(paste0(folder_data, "temp/34-02-contig_size.png"), plot = p, width = 5, height = 5)



























