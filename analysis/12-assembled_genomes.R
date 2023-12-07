#' This script summarize assembled genome information

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(seqinr)
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    mutate(genome_id = factor(genome_id, genome_id))

# 1. aggregate the contig information in the assembled genomes

list_g_contigs <- rep(list(NA), length(isolates$genome_name))
for (i in 1:length(list_g_contigs)) {
    fa <- read.fasta(paste0(folder_genomes, isolates$genome_name[i], "/02-denovo_assembly/genome.fasta"))
    fa_len <- sapply(fa, length)
    list_g_contigs[[i]] <- tibble(genome_id = isolates$genome_id[i], contig_id = names(fa_len), contig_length = fa_len)
}

contigs <- bind_rows(list_g_contigs) %>%
    arrange(genome_id, desc(contig_length)) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))

write_csv(contigs, paste0(folder_data, "temp/12-contigs.csv"))

# Remove small contigs
contigs_large <- contigs %>%
    group_by(genome_id) %>%
    filter(contig_length > 500000)

write_csv(contigs_large, paste0(folder_data, "temp/12-contigs_large.csv"))

