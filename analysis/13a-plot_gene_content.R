#' This script plots the heatmap of gene content

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

genomes <- read_csv(paste0(folder_data, "temp/00-genomes.csv")) 
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
pa <- read_delim(paste0(folder_data, "genomics/pangenome/panaroo/gene_presence_absence.Rtab"))
pa <- pa %>% clean_names()


# 
pal <- pa %>%
    pivot_longer(cols = -gene, names_to = "genome_id") %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    mutate(value = factor(value))

pal_n <- pal %>%
    group_by(gene) %>%
    filter(value == 1) %>%
    count() %>%
    ungroup %>%
    arrange(desc(n))

#
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

ggsave(paste0(folder_data, "temp/13a-01-gene_content.png"), plot = p, width = 15, height = 6)

nrow(pal_n) # 31964 genes in the pangenome
