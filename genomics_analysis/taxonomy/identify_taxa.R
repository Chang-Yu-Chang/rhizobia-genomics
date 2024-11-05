#' This script assigns the taxonomy

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, 'mapping/isolates.csv'))
genomes <- read_csv(paste0(folder_data, 'genomics_analysis/genomes/genomes.csv'))

# 1. Genome blast. Clean the replicon name ----
b_genome <- read_csv(paste0(folder_data, 'genomics_analysis/taxonomy/b_genome.csv'))
contigs <- b_genome %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    mutate(contig_id = paste0(genome_id, "_", qseqid)) %>%
    mutate(species = paste0("E. ", species)) %>%
    mutate(replicon = case_when(
        str_detect(replicon, "chromosome") ~ "chromosome",
        str_detect(replicon, "psymA") ~ str_replace(replicon, "psymA", "pSymA"),
        str_detect(replicon, "psymB") ~ str_replace(replicon, "psymB", "pSymB"),
        str_detect(replicon, "pA") ~ str_replace(replicon, "pA", "pSymA"),
        str_detect(replicon, "pB") ~ str_replace(replicon, "pB", "pSymB"),
        str_detect(replicon, "pSMED01") ~ str_replace(replicon, "pSMED01", "pSymB like"),
        str_detect(replicon, "pSMED02") ~ str_replace(replicon, "pSMED02", "pSymA like"),
        str_detect(replicon, "pSMED03") ~ str_replace(replicon, "pSMED03", "accessory"),
        str_detect(replicon, "accessoryA") ~ str_replace(replicon, "accessoryA", "accessory"),
        str_detect(replicon, "pWSM1115_1") ~ str_replace(replicon, "pWSM1115_1", "pSymB like"),
        str_detect(replicon, "pWSM1115_2") ~ str_replace(replicon, "pWSM1115_2", "pSymA like"),
        str_detect(replicon, "pWSM1115_3") ~ str_replace(replicon, "pWSM1115_3", "accessory"),
        str_detect(replicon, "pSU277_1") ~ str_replace(replicon, "pSU277_1", "pSymB like"),
        str_detect(replicon, "pSU277_2") ~ str_replace(replicon, "pSU277_2", "pSymA like"),
        str_detect(replicon, "pSU277_3") ~ str_replace(replicon, "pSU277_3", "accessory"),
        str_detect(replicon, "pSU277_4") ~ str_replace(replicon, "pSU277_4", "accessory"),
        T ~ replicon
    )) %>%
    # Remove USDA1022 reference with bad assembled contigs
    filter(strain != "USDA1022")

write_csv(contigs, paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))

# 2. 16S blast ----
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
contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))
rrnas <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/rrnas.csv"))

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

colnames(contigs)[2:5] <- paste0("contig_", colnames(contigs)[2:5])
colnames(rrnas)[2:5] <- paste0("rrna_", colnames(rrnas)[2:5])
isolates_tax <- left_join(contigs, rrnas)

write_csv(isolates_tax, paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

# Clean create a master isolate table ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

iso <- isolates %>%
    left_join(isolates_tax) %>%
    mutate(population = ifelse(population == "VA", "elevation", "urbanization")) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id) %>%
    select(population, site_group, exp_id, genome_id, starts_with("rrna"), starts_with("contig")) %>%
    rename(site = site_group, genome = genome_id) %>%
    mutate(contig_length = round(contig_length/1000, 2),
           contig_pident = round(contig_pident, 1),
           rrna_pident = round(rrna_pident, 1)) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything()) %>%
    select(` `, population, site, exp_id, genome, rrna_species, rrna_pident, contig_species, contig_pident)

write_csv(iso, paste0(folder_data, "output/iso.csv"))
