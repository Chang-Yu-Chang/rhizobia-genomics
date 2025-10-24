
library(tidyverse)
source(here::here("metadata.R"))

#blast <- read_tsv(paste0(folder_genomics, "ibd/blast_results/g10_vs_refs.tsv"), col_names = c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore"))
folder_genomics <- "~/Dropbox/lab/rhizobia-genomics/data/genomics"
folder_ibd <- file.path(folder_genomics, "ibd/")
folder_blast <- file.path(folder_ibd, "blast_results/")

# Reference map
replicon_map <- tribble(
    ~sseqid,         ~species,               ~replicon,
    "NC_003047.1",   "S. meliloti",          "chromosome",
    "NC_003037.1",   "S. meliloti",          "pSymA",
    "NC_003078.1",   "S. meliloti",          "pSymB",
    "NC_009636.1",   "S. medicae",           "chromosome",
    "NC_009620.1",   "S. medicae",           "pSymA",
    "NC_009621.1",   "S. medicae",           "pSymB",
    "NC_009622.1",   "S. medicae",           "pAcce"
)

# Blast results
blast_files <- list.files(folder_blast, pattern = "_vs_refs.tsv$", full.names = TRUE)

blast_all <- map_dfr(blast_files, function(f) {
    sample_name <- str_remove(basename(f), "_vs_refs.tsv$")
    read_tsv(f, show_col_types = FALSE,
             col_names = c("qseqid", "sseqid", "qlen", "pident", "length", "evalue", "bitscore")) %>%
        mutate(genome_id = sample_name)
})

# Filter for the top hit
contigs <- blast_all %>%
    left_join(replicon_map, by="sseqid") %>%
    group_by(genome_id, qseqid) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    filter(replicon != "pAcce") %>%
    ungroup()

write_csv(contigs, paste0(folder_ibd, "contigs.csv"))
