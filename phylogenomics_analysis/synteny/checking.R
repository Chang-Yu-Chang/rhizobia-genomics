#' Plots the synteny

renv::load()
library(tidyverse)
library(gggenomes)
source(here::here("metadata.R"))

contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

i=16
isolates$genome_id[i]
m_genes <- read_gff3(paste0(folder_genomics, "gff/genomes/", isolates$genome_id[i], ".gff"))
table(m_genes$type)

string_interest = ""
m_genes %>%
    filter(seq_id == "contig_2") %>%
    filter(!str_detect(product, "hypothetical")) %>%
    view
    filter(str_detect(tolower(product), string_interest) | str_detect(tolower(name), string_interest)) %>%
    select(seq_id, start, end, type, name, product)

isolates_tax %>%
    filter(genome_id == isolates$genome_id[i]) %>%
    select(genome_id, sm_species, contig_species, rrna_species)

contigs_name = "contig_2"
contigs %>%
    filter(genome_id == isolates$genome_id[i]) %>%
    filter(contig_id == paste0(isolates$genome_id[i], "_", contigs_name)) %>%
    select(genome_id, contig_id, species, strain, replicon, pident)
