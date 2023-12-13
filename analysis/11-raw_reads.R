#' This script summarize assembled genome information

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(seqinr)
    source(here::here("analysis/00-metadata.R"))
})


isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    filter(!is.na(rhizobia_population))

# 0. data wrangling ----
# 0.1 join data ----
compute_q <- function (asc) {
    #' Compute the mean phred (Quality score) of a raw read
    #' equation is from https://labs.epi2me.io/quality-scores/
    phred <- as.numeric(charToRaw(asc))-33
    round(-10*log10(sum(10^(-phred/10))/nchar(asc)), 2)
}

list_filtered_reads <- rep(list(NA), nrow(isolates))

for (i in 1:nrow(isolates)) {
    list_filtered_reads[[i]] <- read_table(paste0(folder_genomes, "/", isolates$genome_name[i], "/01-reads_qc/filtered_reads.txt"), col_names = c("name", "asc", "length"), show_col_types = F) %>%
        rowwise() %>%
        mutate(q_score = compute_q(asc), .keep = "unused") %>%
        mutate(genome_id = isolates$genome_id[i]) %>%
        ungroup()
    cat(i)
}


filtered_reads <- bind_rows(list_filtered_reads)

write_csv(filtered_reads, paste0(folder_data, "temp/11-filtered_reads.csv"))

