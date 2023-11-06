#' This script is to have a overview on the raw read statistics

library(tidyverse)
library(cowplot)
library(seqinr)
source(here::here("analysis/00-metadata.R"))

# 0. read the raw read data
if (FALSE) {
    #' Aggregate the contig data
    #' This only needs to run once
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


}
g_contigs <- read_csv(paste0(folder_data, "temp/34-g_contigs.csv"), show_col_types = F)
isolates_rhizo <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F)

# Clean up
g_contigs <- g_contigs %>%
    left_join(isolates_rhizo) %>%
    mutate(genome_id = factor(genome_id, paste0("g", 1:19))) %>%
    filter(genus == "Ensifer")

g_contigs %>%
    group_by(genome_id) %>%
    summarize(count = n())

# 1. plot the genome size
p <- g_contigs %>%
    group_by(genome_id) %>%
    mutate(contig_ordered = factor(1:n())) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = contig_length/10^6, group = contig_ordered, fill = contig_ordered), position = position_stack(reverse = T), color = "black", width = 0.7) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 8.5), breaks = 1:8) +
    scale_fill_manual(values = c("grey", rep("white", 20)) %>% setNames(1:21), breaks = c(1,2), labels = c("largest contig", "other contig")) +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank()
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "genome", y = "genome size (Mbp)")

ggsave(paste0(folder_data, "temp/34-01-genome_size.png"), plot = p, width = 8, height = 4)

# What's the genome size?
g_genomes <- g_contigs %>%
    group_by(genome_id) %>%
    summarize(genome_size = sum(contig_length)/10^6)
g_genomes
# genome_id genome_size
# <fct>           <dbl>
# 1 g2               6.79
# 2 g3               7.51
# 3 g4               7.60
# 4 g5               7.56
# 5 g6               7.59
# 6 g8               7.52
# 7 g9               7.01
# 8 g10              7.58
# 9 g11              7.31
# 10 g13              7.43
# 11 g15              8.03
# 12 g16              7.41
# 13 g17              7.42
# 14 g19              7.18

round(range(g_genomes$genome_size), 2) # [6.79, 8.03] Mbp

# 2. plot the contig size for all genomes
p <- g_contigs %>%
    ggplot() +
    geom_histogram(aes(x = contig_length/10^6), color = "black", fill = "white") +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "contig size (Mbp)", y = "count")

ggsave(paste0(folder_data, "temp/34-02-contig_size.png"), plot = p, width = 4, height = 4)



























