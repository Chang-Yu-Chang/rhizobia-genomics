#' This script plots the gene content by the contigs

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    library(ggsci)
    library(proxy) # for computing jaccard matrix
    library(ape)
    library(ggtree)
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    filter(!genome_id %in% paste0("g", c(1,7,12,14,18)))
egc <- read_csv(paste0(folder_data, "temp/17-egc.csv"), show_col_types = F) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))
gcalls <- read_csv(paste0(folder_data, "temp/17-gcalls.csv"), show_col_types = F)
egcalls <- read_csv(paste0(folder_data, "temp/17-egcalls.csv"), show_col_types = F)
nrow(isolates) # This should be 32 ensifer isolates + 4 ncbi genomes


# 1. plot the copy number variation ----
egc_copy <- egc %>%
    group_by(genome_id, gene_cluster_id) %>%
    count(name = "n_copies")

table(egc_copy$n_copies)
p <- egc_copy %>%
    #filter(n_copies <= 10) %>%
    #filter(genome_id == "em1021") %>%
    ggplot() +
    geom_histogram(aes(x = n_copies), color = "black", fill = "white", binwidth = 1) +
    scale_x_continuous(breaks = c(1:10, seq(5, 50, 5))) +
    scale_y_log10() +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/17c-01-copy_number.png"), p, width = 5, height = 5)

# 2. plot the copy number variation in each sample ----
p <- egc_copy %>%
    #filter(n_copies <= 10) %>%
    #filter(genome_id == "em1021") %>%
    ggplot() +
    geom_histogram(aes(x = n_copies), color = "black", fill = "white", binwidth = 1) +
    scale_y_log10() +
    facet_wrap(~genome_id, ncol = 6) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/17c-02-copy_number.png"), p, width = 10, height = 10)


# 3. find the gene with the highest copy number ----
egc_copy %>%
    filter(n_copies > 10) %>%
    ungroup() %>%
    distinct(gene_cluster_id)

#














