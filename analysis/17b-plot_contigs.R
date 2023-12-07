#' This script plots the gene content by the contigs

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    library(ggsci)
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    filter(!genome_id %in% paste0("g", c(1,7,12,14,18)))
egc <- read_csv(paste0(folder_data, "temp/17-egc.csv"), show_col_types = F)
gcalls <- read_csv(paste0(folder_data, "temp/17-gcalls.csv"), show_col_types = F)
contigs_mash <- read_csv(paste0(folder_data, "temp/14-contigs_mash.csv"), show_col_types = F)
nrow(isolates) # This should be 32 ensifer isolates + 4 ncbi genomes

#n_g <- nrow(isolates)

# 0. join data ----
egcalls <- egc %>%
    left_join(gcalls) %>%
    select(unique_id, gene_cluster_id, genome_id, contig_id, gene_callers_id, max_num_paralogs, prokka_prodigal_acc) %>%
    left_join(contigs_mash) %>%
    right_join(isolates)


# 1. plot the number of genes on each contig ----
p <- egcalls %>%
    distinct(genome_id, contig_id, contig_length, prokka_prodigal_acc) %>%
    group_by(genome_id, contig_id, contig_length) %>%
    count(name = "n_genes") %>%
    mutate(contig_length = contig_length/10^6,
           n_genes = n_genes/10^3) %>%
    ggplot() +
    geom_point(aes(x = contig_length, y = n_genes), shape = 21) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "contig length (Mbp)", y = "# of genes (k)")
ggsave(paste0(folder_data, "temp/17b-01-contig_length_vs_genes.png"), p, width = 4, height = 4)

# 2. plot the
