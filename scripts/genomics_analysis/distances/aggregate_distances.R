#' This script

library(tidyverse)
library(janitor)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

# ANI
dist_ani <- read_delim(paste0(folder_data, "genomics/distances/ani/ani_genomes.txt"), delim = "\t", col_names = F)
colnames(dist_ani) <- c("genome_id1", "genome_id2", "d_ani", "frag1", "frag2")

dist_ani <- dist_ani %>%
    # Compute so that d_ani represents the genomic distance instance of similarity
    mutate(d_ani = 1-d_ani/100) %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = ordered(genome_id1, isolates$genome_id), genome_id2 = ordered(genome_id2, isolates$genome_id)) %>%
    arrange(genome_id1, genome_id2) %>%
    select(genome_id1, genome_id2, d_ani)

# k-mers
dist_kmer <- read_delim(paste0(folder_data, "genomics/distances/kmer/kmer.txt"))
list_sigs <- read_delim(paste0(folder_data, "genomics/distances/kmer/list_sigs.txt"), delim = "\t", col_names = F)

dist_kmer <- dist_kmer %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_kmer")  %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = ordered(genome_id1, isolates$genome_id), genome_id2 = ordered(genome_id2, isolates$genome_id)) %>%
    arrange(genome_id1, genome_id2)


# Join the long format distance tables
distl <- dist_ani %>%
    left_join(dist_kmer)

# Save the files
nrow(distl) # 38^2 = 1444
write_csv(distl, paste0(folder_data, "genomics_analysis/distances/distl.csv"))
