#' This script joins contig data

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))
top_genome <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/top_genome.csv"))


# Contig level taxonomy
# b_genome <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/b_genome.csv")) %>%
#     rename(contig_id = qseqid) %>%
#     mutate(contig_id = paste0(genome_id, "_", contig_id)) %>%
#     select(genome_id, contig_id, species, strain, replicon, pident, length)
#
# blast_genome <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/blast_genome.csv")) %>%
#     rename(contig_id = qseqid) %>%
#     mutate(contig_id = paste0(genome_id, "_", contig_id))
#     select(genome_id, contig_id, species, strain, replicon, pident, length) %>%
#     # Remove those that are symbiotic species but blast to nonsymbiotic plasmids
#     filter(!(!genome_id %in% c("g2", "g3", "g15") & species %in% c("canadensis", "adhaerens")))
#     top_genome <- blast_genome %>%
#     group_by(genome_id, contig_id) %>%
#     arrange(desc(bitscore)) %>%
#     # Find the top hit
#     slice(1)


# Plasmid
plasmids <- tibble(
    species = c(rep("meliloti", 3*4+1), rep("medicae", 3*4+1), "canadensis", "adhaerens"),
    strain = c(rep(c("USDA1106", "1021", "WSM1022"), each = 3), rep(c("USDA1021", "WSM419", "WSM1115", "SU277"), each = 4), "SU277", "T173", "Corn53"),
    replicon_type = c(rep(c("chromosome", "pSymA", "pSymB"), 3), rep(c("chromosome", "pSymA", "pSymB", "pAcce"), 4), "pAcce", "chromosome", "chromosome"),
    replicon = c(
        "chromosome", "psymA", "psymB",
        "chromosome", "pSymA", "pSymB",
        "chromosome", "pA", "pB",
        "chromosome", "psymA", "psymB", "accessoryA",
        "chromosome", "pSMED02", "pSMED01", "pSMED03",
        "chromosome", "pWSM1115_2", "pWSM1115_1", "pWSM1115_3",
        "chromosome", "pSU277_2", "pSU277_1", "pSU277_3", "pSU277_4",
        "chromosome", "chromosome"
    )
)

# Join the datasets
contigs <- genomes %>%
    left_join(top_genome) %>%
    mutate(replicon = str_remove(replicon, "plasmid ")) %>%
    mutate(frac_blas = round(length / contig_length * 100, 2))

#
contigs <- contigs %>%
    drop_na(species) %>%
    filter(contig_length > 1e4) %>%
    left_join(plasmids)


write_csv(contigs, paste0(folder_data, "genomics_analysis/genomes/contigs.csv"))
