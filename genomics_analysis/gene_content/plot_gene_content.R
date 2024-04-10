#' This script plots the heatmap of gene content

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))


#
pal <- gpa %>%
    pivot_longer(cols = -gene, names_to = "genome_id") %>%
    filter(genome_id %in% isolates$genome_id) %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    mutate(value = factor(value))

paln <- pal %>%
    group_by(gene) %>%
    filter(value == 1) %>%
    count() %>%
    ungroup %>%
    arrange(desc(n))

# 1. Heatmap
p <- pal %>%
    mutate(gene = factor(gene, pal_n$gene)) %>%
    ggplot() +
    geom_tile(aes(x = gene, y = genome_id, fill = value)) +
    scale_fill_manual(values = c(`1` = "maroon", `0` = "white")) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    guides(fill = "none") +
    labs()

ggsave(paste0(folder_data, "genomics_analysis/gene_content/01-heatmap.png"), p, width = 15, height = 6)

nrow(pal_n) # 31964 genes in the pangenome

# 2. gene frequency spectrum
p <- paln %>%
    ggplot() +
    geom_histogram(aes(x = n), color = "black", fill = "white") +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "# of genomes")

ggsave(paste0(folder_data, "genomics_analysis/gene_content/02-gene_frequency_spectrum.png"), p, width = 5, height = 4)






















