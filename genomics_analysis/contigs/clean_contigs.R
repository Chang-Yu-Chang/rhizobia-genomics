#' This script joins contig data

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

# Replicon size by nt
genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))

# Contig level taxonomy
b_genome <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/b_genome.csv")) %>%
    rename(contig_id = qseqid) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id)) %>%
    select(genome_id, contig_id, species, strain, replicon, pident, length)

# Plasmid
plasmids <- tibble(
    species = c(rep("meliloti", 4*3), rep("medicae", 3*3), "canadensis", "adhaerens"),
    strain = c(rep(c("USDA1106", "1021", "WSM1022", "USDA1021", "WSM419", "WSM1115", "SU277"), each = 3), "T173", "Corn53"),
    replicon_type = c(rep(c("chromosome", "psymA like", "psymB like"), 7), "chromosome", "chromosome"),
    replicon = c(
        "chromosome", "psymA", "psymB",
        "chromosome", "pSymA", "pSymB",
        "chromosome", "pA", "pB",
        "chromosome", "psymA", "psymB",
        "chromosome", "pSMED02", "pSMED01",
        "chromosome", "pWSM1115_2", "pWSM1115_1",
        "chromosome", "pSU277_2", "pSU277_1",
        "chromosome", "chromosome"
    )
)

# Join the datasets
contigs <- genomes %>%
    left_join(b_genome) %>%
    mutate(replicon = str_remove(replicon, "plasmid ")) %>%
    mutate(frac_blas = round(length / contig_length * 100, 2))

#
contigs <- contigs %>%
    drop_na(species) %>%
    filter(contig_length > 1e6) %>%
    left_join(plasmids) %>%
    #drop_na(replicon_type) %>%
    # remove odd contigs
    filter(!contig_id %in% c("g28_contig_15"))
# g20_contig_3  and g20_contig7 together is chromosome


write_csv(contigs, paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))
