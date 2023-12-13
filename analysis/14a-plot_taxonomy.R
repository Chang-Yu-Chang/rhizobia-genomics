#' This script assigns the taxonomy to isolates

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F)
contigs_large <- read_csv(paste0(folder_data, "temp/12-contigs_large.csv"), show_col_types = F)
contigs <- read_csv(paste0(folder_data, "temp/12-contigs.csv"), show_col_types = F)
blast_genomes_top <- read_csv(paste0(folder_data, "temp/14-blast_genomes_top.csv"), show_col_types = F)
blast_16s <- read_csv(paste0(folder_data, "temp/14-blast_16s.csv"), show_col_types = F)


# 0. join the data ----
contigs_blast <- blast_genomes_top %>%
    rename(contig_id = qseqid) %>%
    left_join(contigs_large) %>%
    filter(genome_id %in% isolates$genome_id[!is.na(isolates$exp_id)]) %>%
    drop_na


# 1. plot the tree based on genome ----
