#' This script assigns the taxonomy

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, 'mapping/isolates.csv'))
genomes <- read_csv(paste0(folder_data, 'genomics_analysis/genomes/genomes.csv'))

# 1. sourmash ----
sm_genome <- read_csv(paste0(folder_data, 'genomics_analysis/taxonomy/sm_genome.csv'))
sms <- sm_genome %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    group_by(genome_id) %>%
    mutate(species = str_split(name, " ") %>% `[[`(1) %>% `[`(c(2,3)) %>% paste(collapse = " ")) %>%
    mutate(species = str_replace(species, "Sinorhizobium", "Ensifer")) %>%
    select(genome_id, species, query_containment_ani, name)
write_csv(sms, paste0(folder_data, "genomics_analysis/taxonomy/sms.csv"))

# 2. Genome blast. Clean the replicon name ----
b_genome <- read_csv(paste0(folder_data, 'genomics_analysis/taxonomy/b_genome.csv'))
contigs <- b_genome %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    mutate(contig_id = paste0(genome_id, "_", qseqid)) %>%
    mutate(species = paste0("E. ", species)) %>%
    #select(genome_id, contig_id, species, strain, replicon, bitscore) %>%
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

# 3. 16S blast ----
b_16s <- read_csv(paste0(folder_data, 'genomics_analysis/taxonomy/b_16s.csv'))
rrnas <- b_16s %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    group_by(genome_id) %>%
    mutate(species = paste0("E. ", str_split(scomment, " ")[[1]] %>% `[`(2))) %>%
    select(genome_id, species, sseqid, pident, length) %>%
    ungroup()
write_csv(rrnas, paste0(folder_data, "genomics_analysis/taxonomy/rrnas.csv"))

# Assign isolates to taxonomy
genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))
sms <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/sms.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))
rrnas <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/rrnas.csv"))

sms <- sms %>%
    group_by(genome_id) %>%
    arrange(desc(query_containment_ani)) %>%
    slice(1) %>%
    select(-name)
contigs <- contigs %>%
    group_by(genome_id) %>%
    left_join(genomes) %>%
    filter(replicon == "chromosome") %>%
    arrange(genome_id, desc(contig_length)) %>%
    slice(1) %>%
    select(genome_id, species, sseqid, pident, length) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id) %>%
    ungroup()
rrnas <- rrnas %>%
    group_by(genome_id) %>%
    arrange(desc(pident)) %>%
    slice(1) %>%
    ungroup()

colnames(sms)[2:3] <- paste0("sm_", colnames(sms)[2:3])
colnames(contigs)[2:5] <- paste0("contig_", colnames(contigs)[2:5])
colnames(rrnas)[2:5] <- paste0("rrna_", colnames(rrnas)[2:5])

isolates_tax <- rrnas %>%
    left_join(contigs)

write_csv(isolates_tax, paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
