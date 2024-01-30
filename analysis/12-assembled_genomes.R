#' This script aggregates the assembled genome fasta

renv::load()
library(tidyverse)
library(janitor)
library(seqinr)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv")) 
isolates <- isolates %>% drop_na(exp_id)

# Aggregate the contig information in the assembled genomes
list_g_contigs <- rep(list(NA), nrow(isolates))
for (i in 1:length(list_g_contigs)) {
    fa <- read.fasta(paste0(folder_genomics, "/genomes/", isolates$genome_id[i], ".fasta"))
    fa_len <- sapply(fa, length)
    list_g_contigs[[i]] <- tibble(genome_id = isolates$genome_id[i], contig_id = names(fa_len), contig_length = fa_len)
    cat(i)
}

contigs <- bind_rows(list_g_contigs[!is.na(list_g_contigs)]) %>%
    # Arrange contig by length
    arrange(genome_id, desc(contig_length)) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    # Remove small contigs < 10kb
    group_by(genome_id) %>%
    filter(contig_length > 10000)

write_csv(contigs, paste0(folder_data, "temp/12-contigs.csv"))

