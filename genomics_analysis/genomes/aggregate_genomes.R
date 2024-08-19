#' This script aggregates the assembled genome fasta

renv::load()
library(tidyverse)
library(janitor)
library(seqinr)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

# Aggregate the contig information in the assembled genomes ----
list_g_contigs <- rep(list(NA), nrow(isolates))
for (i in 1:length(list_g_contigs)) {
    # This is reading the fasta in the consolidated genome folder
    fa <- read.fasta(paste0(folder_genomics, "fasta/genomes/", isolates$genome_id[i], ".fasta"))
    fa_len <- sapply(fa, length)
    list_g_contigs[[i]] <- tibble(genome_id = isolates$genome_id[i], contig_id = names(fa_len), contig_length = fa_len)
    cat("\t", i)
}

genomes <- bind_rows(list_g_contigs) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id)) %>%
    # Arrange contigs by length
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id, desc(contig_length)) %>%
    # Remove small contigs < 100kb
    filter(contig_length > 100000)

write_csv(genomes, paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))

# Extract QC values ----
## quast
list_g <- rep(list(NA), nrow(isolates))
for (i in 1:length(list_g)) {
    rp <- read_delim(paste0(folder_genomics, "assembly/", isolates$genome_id[i],"/quast/report.tsv")) %>%
        set_names(c("name", "value"))

    list_g[[i]] <- tibble(
        genome_id = isolates$genome_id[i],
        name = c("GC (%)", "Largest contig", "Total length (>= 10000 bp)", "# contigs (>= 10000 bp)", "N50", "N90", "L50", "L90")
    ) %>%
        left_join(rp)
}

quast <- bind_rows(list_g) %>%
    group_by(genome_id) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    clean_names()

# ## busco
# list_b <- rep(list(NA), nrow(isolates))
# for (i in 1:length(list_b)) {
#     file_name <- paste0(folder_genomics, "assembly/", isolates$genome_id[i],"/busco/run_alphaproteobacteria_odb10/short_summary.txt")
#     rls <- readLines(file_name)
#     tmp <- rls[10:15] %>% str_split("\t", simplify = T) %>% `[`(,2) %>% as.numeric()
#     list_b[[i]] <- tibble(genome_id = isolates$genome_id[i], name = c("C", "S", "D", "F", "M", "T"), value = tmp)
# }
#
# busco <- bind_rows(list_b) %>%
#     group_by(genome_id) %>%
#     pivot_wider(names_from = name, values_from = value) %>%
#     mutate(completeness = C/T)
#
# qcs <- left_join(quast, busco)
qcs <- quast
write_csv(qcs, paste0(folder_data, "genomics_analysis/genomes/qcs.csv"))







