#' This script clean the gene content data

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv")) 
pa <- read_delim(paste0(folder_data, "genomics/pangenome/panaroo/gene_presence_absence.Rtab"))
pa <- pa %>% clean_names()
dist_fluidity <- read_csv(paste0(folder_data, "temp/13-dist_fluidity.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))


# Remove two genomes with bad assembly quality
isolates_sp <- isolates_contigs %>% 
    filter(!genome_id %in% c("g20", "g28")) %>%
    distinct(genome_id, species) 

# Calculate the genomic fluidifty of all strains
list_ensifer <- isolates_sp$genome_id
length(list_ensifer) # 40 ensifer strain
f_all <- dist_fluidity %>%
    filter(genome_id1 %in% list_ensifer, genome_id2 %in% list_ensifer) %>%
    # Choose only upper triangle
    mutate(genome_id1 = ordered(genome_id1, isolates$genome_id), genome_id2 = ordered(genome_id2, isolates$genome_id)) %>%
    filter(genome_id1 < genome_id2) 
nrow(f_all) # choose(40,2) = 780
sum(f_all$d_fluidity) / choose(length(list_ensifer),2) # 0.354


# Calculate the genomic fluidity of meliloti
list_meliloti <- isolates_sp$genome_id[isolates_sp$species == "meliloti"]
length(list_meliloti) # 25 meliloti
f_meliloti <- dist_fluidity %>%
    filter(genome_id1 %in% list_meliloti & genome_id2 %in% list_meliloti) %>%
    # Choose only upper triangle
    mutate(genome_id1 = ordered(genome_id1, isolates$genome_id), genome_id2 = ordered(genome_id2, isolates$genome_id)) %>%
    filter(genome_id1 < genome_id2) 
nrow(f_meliloti) # choose(25,2) = 300
sum(f_meliloti$d_fluidity) / choose(length(list_meliloti),2) # 0.29


# Calculate the genomic fluidity of medicae
list_medicae <- isolates_sp$genome_id[isolates_sp$species == "medicae"]
length(list_medicae) # 15 meliloti
f_medicae <- dist_fluidity %>%
    filter(genome_id1 %in% list_medicae & genome_id2 %in% list_medicae) %>%
    # Choose only upper triangle
    mutate(genome_id1 = ordered(genome_id1, isolates$genome_id), genome_id2 = ordered(genome_id2, isolates$genome_id)) %>%
    filter(genome_id1 < genome_id2) 
nrow(f_medicae) # choose(15,2) = 105
sum(f_medicae$d_fluidity) / choose(length(list_medicae),2) # 0.187