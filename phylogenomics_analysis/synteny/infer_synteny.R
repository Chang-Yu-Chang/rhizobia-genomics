#' This script inferes the synteny

renv::load()
library(tidyverse)
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))

#
contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))
contigs <- contigs %>%
    select(genome_id, contig_id, species, replicon)


#
db <- read_csv(paste0(folder_genomics, "pangenome/isolates/gene_data.csv"))
db <- db %>%
    rename(genome_id = gff_file, contig_id = scaffold_name) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id))

db %>%
    group_by(genome_id, contig_id) %>%
    count() %>%
    left_join(genomes) %>%
    left_join(contigs) %>%
    filter(contig_length > 10^6) %>%
    view


