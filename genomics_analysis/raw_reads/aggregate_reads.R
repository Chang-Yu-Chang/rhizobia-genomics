#' This script aggregates filtered reads

renv::load()
library(tidyverse)
library(janitor)
library(seqinr)
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))

# Join data
compute_q <- function (asc) {
    #' Compute the mean phred (Quality score) of a raw read
    #' equation is from https://labs.epi2me.io/quality-scores/
    phred <- as.numeric(charToRaw(asc))-33
    round(-10*log10(sum(10^(-phred/10))/nchar(asc)), 2)
}

list_filtered_reads <- rep(list(NA), 32)

for (i in 1:32) {
    list_filtered_reads[[i]] <- read_table(paste0(folder_genomics, "/assembly/", genomes$genome_id[i], "/filtered_reads.txt"), col_names = c("name", "asc", "length"), show_col_types = F) %>%
        rowwise() %>%
        mutate(q_score = compute_q(asc), .keep = "unused") %>%
        mutate(genome_id = genomes$genome_id[i]) %>%
        ungroup()
    cat(i)
}

filtered_reads <- bind_rows(list_filtered_reads)

write_csv(filtered_reads, paste0(folder_data, "genomics_analysis/raw_reads/filtered_reads.csv"))
