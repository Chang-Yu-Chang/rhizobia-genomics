#' This script aggregates the genomes

renv::load()
library(tidyverse)
library(gggenomes)
source(here::here("metadata.R"))


# Aggregate files
contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/contigs.csv"))
detect_replicon <- function (contigs, rr = "pSymA") {
    contigs %>%
        filter(genome_id == isolates$genome_id[i]) %>%
        filter(str_detect(replicon, rr)) %>%
        pull(qseqid)
}
filter_seqs <- function (seqs, length_min = 1e6, rr_con) {
    seqs %>%
        filter(seq_id %in% rr_con) %>%
        filter(length > length_min) %>%
        mutate(seq_id = paste0(isolates$genome_id[i], "_", seq_id))
}
filter_genes <- function (genes, rr_con) {
    m_genes %>%
        filter(seq_id %in% rr_con) %>%
        #filter(str_detect(name, regex)) %>%
        mutate(seq_id = paste0(isolates$genome_id[i], "_", seq_id))
}

tb <- tibble(genome_id = isolates$genome_id, seqs_chr = NA, genes_chr = NA, seqs_a = NA, genes_a = NA, seqs_b = NA, genes_b = NA)
for (i in 1:nrow(tb)) {
#for (i in c(2,3)) {
    if (i == 23) next
    #
    m_seqs <- read_seqs(paste0(folder_genomics, "fasta/genomes/", isolates$genome_id[i], ".fasta"))
    m_genes <- read_gff3(paste0(folder_genomics, "gff/genomes/", isolates$genome_id[i], ".gff"))

    # psyma
    psyma_con <- detect_replicon(contigs, "pSymA")
    tb$seqs_a[i] <- list(filter_seqs(m_seqs, 1e6, psyma_con))
    tb$genes_a[i] <- list(filter_genes(m_genes, psyma_con))

    # psymb
    psymb_con <- detect_replicon(contigs, "pSymB")
    tb$seqs_b[i] <- list(filter_seqs(m_seqs, 1e6, psymb_con))
    tb$genes_b[i] <- list(filter_genes(m_genes, psymb_con))

    # chromosome
    chrom_con <- detect_replicon(contigs, "chromosome")
    tb$seqs_chr[i] <- list(filter_seqs(m_seqs, 1e6, chrom_con))
    tb$genes_chr[i] <- list(filter_genes(m_genes, chrom_con))

    # tb$seqs[i] <- list(ms)
    # tb$genes[i] <- list(mg)
}

tb <- filter(tb, genome_id != "g28")
#tb <- tb[c(2,3),]
mss <- bind_rows(
    mutate(bind_rows(tb$seqs_chr), replicon = "chromosome"),
    mutate(bind_rows(tb$seqs_a), replicon = "psyma"),
    mutate(bind_rows(tb$seqs_b), replicon = "psymb")
)

mgs <- bind_rows(
    mutate(bind_rows(tb$genes_chr), replicon = "chromosome"),
    mutate(bind_rows(tb$genes_a), replicon = "psyma"),
    mutate(bind_rows(tb$genes_b), replicon = "psymb")
)

save(mss, mgs, file = paste0(folder_data, "phylogenomics_analysis/synteny/gggenomes.rdata"))
# write_csv(mss, paste0(folder_data, "phylogenomics_analysis/synteny/mss.csv"))
# write_csv(mgs, paste0(folder_data, "phylogenomics_analysis/synteny/mgs.csv"))









