#' This script plots the contigs

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
contigs <- read_csv(paste0(folder_data, "temp/12-contigs.csv")) %>%
    left_join(genomes) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))


# 1. Number of contigs
p <- contigs %>%
    group_by(genome_id) %>%
    dplyr::count(name = "n_contigs") %>%
    ggplot() +
    geom_histogram(aes(x = n_contigs), color = "black", fill = "white", binwidth = 1) +
    scale_x_continuous(breaks = seq(0, 12, 1)) +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "# of contigs", title = "contigs > 10kb")

ggsave(paste0(folder_data, "temp/12a-01-n_contigs.png"), plot = p, width = 4, height = 4)


# 2. plot the genome size by contigs
p <- contigs %>%
    group_by(genome_id) %>%
    mutate(contig_ordered = factor(1:n())) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = contig_length/10^6, group = contig_ordered, fill = contig_ordered), position = position_stack(reverse = T), color = "black", width = 0.7) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 8.5), breaks = 1:8) +
    scale_fill_manual(values = c("grey", rep("white", 20)) %>% setNames(1:21), breaks = c(1,2), labels = c("largest contig", "other contig")) +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "genome", y = "genome size (Mbp)", title = "contigs > 10kb")

ggsave(paste0(folder_data, "temp/12a-02-genome_size.png"), plot = p, width = 15, height = 6)

# contigs
c_size <- contigs %>%
    left_join(isolates) %>%
    filter(population %in% c("VA", "PA")) %>%
    group_by(genome_id) %>%
    mutate(contig_size = contig_length/10^6)

c_size %>%
    filter(contig_length > 500000) %>%
    filter(genome_id != "g28") %>%
    count() %>%
    pull(n) %>%
    range() # 3 5


# genome size
g_size <- c_size %>%
    group_by(genome_id) %>%
    summarise(genome_size = sum(contig_length)/10^6)

nrow(g_size) # 32 genomes
g_size <- g_size %>% filter(genome_size < 14)
round(range(g_size$genome_size), 2) #   6.79 8.02
round(median(g_size$genome_size), 2) # 7.25













