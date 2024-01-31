#' This script assigns the taxonomy 

renv::load()
library(tidyverse)
library(janitor)
library(seqinr)
library(ape)
source(here::here("analysis/00-metadata.R"))

# Read the list of genomes
list_gff <- read_table(paste0(folder_data, "genomics/pangenome/panaroo/list_gffs.txt"), col_names = F)
genome_ids <- list_gff$X1 %>% str_remove(".+/gff/") %>% str_remove(".gff")
length(genome_ids) # 41 genomes

# Read the list of core genes
folder_alignment <- paste0(folder_data, "genomics/pangenome/panaroo/aligned_gene_sequences/")
list_cg <- list.files(folder_alignment)

# Subset the core genes to only single-copy core genes
tb_core <- tibble(gene_name = list_cg, is_single_copy = NA)
for (i in 1:length(list_cg)) {
    aln <- read.alignment(paste0(folder_alignment, list_cg[i]), format = "fasta")
    list_seq_i <- str_extract(aln$nam, "\\w+;") %>% str_remove(";") %>% str_remove("_R_")
    tb_core$is_single_copy[i] <- length(unique(list_seq_i)) == length(list_seq_i)
    print(i)
}
sum(tb_core$is_single_copy) # 825 single-copy core genes
list_sccg <- list_cg[tb_core$is_single_copy]

write_csv(tibble(x = list_sccg), paste0(folder_data, "temp/15-list_sccg.csv"))

# Compute the hamming distance 
hamming <- function (seq1, seq2) {
    d_ham <- sum(strsplit(seq1, '')[[1]] != strsplit(seq2, '')[[1]])
    return(d_ham / nchar(seq1))
}

tb_sccg <- tibble(gene_name = list_sccg, dat = NA)
for (i in 1:length(list_sccg)) {
    aln <- read.alignment(paste0(folder_alignment, list_sccg[i]), format = "fasta")
    genome_ids <- str_extract(aln$nam, "\\w+;") %>% str_remove(";") %>% str_remove("_R_")
    tb_ham <- crossing(seq_id1 = 1:aln$nb, seq_id2 = 1:aln$nb, d_hamming = NA) %>%
        filter(seq_id1 < seq_id2) %>%
        # Compute hamming
        rowwise() %>%
        mutate(d_hamming = hamming(aln$seq[[seq_id1]], aln$seq[[seq_id2]])) %>%
        ungroup() %>%
        # Mathc genome ids
        mutate(genome_id1 = genome_ids[seq_id1], genome_id2 = genome_ids[seq_id2]) %>%
        select(genome_id1, genome_id2, d_hamming)
    tb_sccg$dat[[i]] <- list(tb_ham)
    print(i)
}

dist_sccg_ham <- tb_sccg %>%
    unnest(dat) %>%
    unnest(dat) %>%
    mutate(genome_id1 = ordered(genome_id1, genome_ids), genome_id2 = ordered(genome_id2, genome_ids)) 
    
to_swap <- which(dist_sccg_ham$genome_id1 > dist_sccg_ham$genome_id2) 
x <- dist_sccg_ham$genome_id2[to_swap]
dist_sccg_ham$genome_id2[to_swap] <- dist_sccg_ham$genome_id1[to_swap]
dist_sccg_ham$genome_id1[to_swap] <- x

dist_sccg_ham %>%
    distinct(genome_id1, genome_id2) %>%
    nrow() # choose(41,2) = 820




# Average across single copy core genes
dist_sccg <- dist_sccg_ham %>%
    group_by(genome_id1, genome_id2) %>%
    summarize(d_hamming = mean(d_hamming))

write_csv(dist_sccg, paste0(folder_data, "temp/15-dist_sccg.csv"))
