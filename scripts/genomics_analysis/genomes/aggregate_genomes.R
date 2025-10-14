#' This script aggregates the assembled genome fasta

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
    arrange(genome_id, desc(contig_length))

write_csv(genomes, paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))

# flye ----
list_flye <- rep(list(NA), nrow(isolates))
for (i in 1:length(list_flye)) {
    list_flye[[i]] <- read_table(paste0(folder_genomics, "assembly/", isolates$genome_id[i], "/flye/assembly_info.txt"), show_col_types = F) %>%
        mutate(genome_id = isolates$genome_id[i]) %>%
        filter(length > 10000) %>%
        select(-graph_path)
}
flye <- bind_rows(list_flye) %>% filter(length > 10^6)
flye %>%
    group_by(genome_id) %>%
    summarize(n = n()) # number of genomes with >= 3 large contigs
sum(flye$circ. == "Y") / nrow(flye) # fraction of ciruclar contigs

# quast ----
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

write_csv(quast, paste0(folder_data, "genomics_analysis/genomes/quast.csv"))

# busco ----
list_b <- rep(list(NA), nrow(isolates))
lins <- tibble(taxlevel = c("class", "order", "family", "genus"), lineage = c("alphaproteobacteria_odb10", "rhizobiales_odb10", "rhizobiaceae_odb12", "sinorhizobium_odb12"))

for (i in 1:length(list_b)) {
    list_lin <- list()
    for (lin in lins$lineage) {
        file_name <- paste0(folder_genomics, "assembly/", isolates$genome_id[i],"/busco/", lin, "/run_", lin,"/short_summary.txt")
        rls <- readLines(file_name)
        tmp <- rls[10:15] %>% str_split("\t", simplify = T) %>% `[`(,2) %>% as.numeric()
        list_lin[lin] <- list(tibble(
            lineage = lin,
            genome_id = isolates$genome_id[i],
            name = c("C", "S", "D", "F", "M", "T"),
            value = tmp
        ))
    }
    list_b[[i]] <- bind_rows(list_lin)
}


busco <- bind_rows(list_b) %>%
    left_join(lins) %>%
    group_by(genome_id, taxlevel, lineage) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    mutate(completeness = C/T)

write_csv(busco, paste0(folder_data, "genomics_analysis/genomes/busco.csv"))

# checkm ----
checkm <- read_lines(paste0(folder_genomics, "checkm/checkm_summary.tsv")) %>%
    discard(~ str_detect(.x, "INFO|---|CheckM|^$")) %>%
    # Split on one or more spaces
    str_squish() %>%
    str_remove(" \\(UID\\d+\\)") %>%
    `[`(c(-1, -20)) %>%
    str_split("\\s", simplify = TRUE) %>%
    as_tibble(.name_repair = ~ paste0("V", seq_along(.)))

checkm <- checkm %>%
    select(V1, V2, V6, V7) %>%
    rename(
        bin_id = V1,
        marker_lineage = V2,
        completeness = V6,
        contamination = V7
    )

write_csv(checkm, paste0(folder_data, "genomics_analysis/genomes/checkm.csv"))
