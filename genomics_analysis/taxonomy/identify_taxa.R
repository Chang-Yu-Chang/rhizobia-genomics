#' This script assigns the taxonomy

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, 'genomics_analysis/genomes/genomes.csv'))
b_genome <- read_csv(paste0(folder_data, 'genomics_analysis/taxonomy/b_genome.csv'))

# Clean the replicon name
contigs <- b_genome %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    mutate(contig_id = paste0(genome_id, "_", qseqid)) %>%
    select(genome_id, contig_id, species, strain, replicon, bitscore) %>%
    mutate(replicon = case_when(
        str_detect(replicon, "chromosome") ~ "chromosome",
        str_detect(replicon, "psymA") ~ str_replace(replicon, "psymA", "pSymA"),
        str_detect(replicon, "psymB") ~ str_replace(replicon, "psymB", "pSymB"),
        str_detect(replicon, "pA") ~ str_replace(replicon, "pA", "pSymA"),
        str_detect(replicon, "pB") ~ str_replace(replicon, "pB", "pSymB"),
        str_detect(replicon, "pSMED01") ~ str_replace(replicon, "pSMED01", "pSymA like"),
        str_detect(replicon, "pSMED02") ~ str_replace(replicon, "pSMED02", "pSymB like"),
        str_detect(replicon, "pSMED03") ~ str_replace(replicon, "pSMED03", "accessory"),
        str_detect(replicon, "accessoryA") ~ str_replace(replicon, "accessoryA", "accessory"),
        str_detect(replicon, "pWSM1115_1") ~ str_replace(replicon, "pWSM1115_1", "pSymA like"),
        str_detect(replicon, "pWSM1115_2") ~ str_replace(replicon, "pWSM1115_2", "pSymB like"),
        str_detect(replicon, "pWSM1115_3") ~ str_replace(replicon, "pWSM1115_3", "accessory"),
        str_detect(replicon, "pSU277_1") ~ str_replace(replicon, "pSU277_1", "pSymA like"),
        str_detect(replicon, "pSU277_2") ~ str_replace(replicon, "pSU277_2", "pSymB like"),
        str_detect(replicon, "pSU277_3") ~ str_replace(replicon, "pSU277_3", "accessory"),
        str_detect(replicon, "pSU277_4") ~ str_replace(replicon, "pSU277_4", "accessory"),
        T ~ replicon
    )) %>%
    # Remove USDA1022 reference with bad assembled contigs
    filter(strain != "USDA1022")

write_csv(contigs, paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))

# Assign isolates to taxonomy
genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))
isolates_contigs <- contigs %>%
    filter(replicon == "chromosome") %>%
    left_join(genomes) %>%
    group_by(genome_id) %>%
    arrange(genome_id, desc(contig_length)) %>%
    slice(1) %>%
    select(genome_id, species, strain, contig_id, contig_length) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id) %>%
    ungroup()

write_csv(isolates_contigs, paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
