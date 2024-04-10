#' This script assigns the taxonomy

renv::load()
library(tidyverse)
library(janitor)
library(seqinr)
library(ape)
source(here::here("metadata.R"))

# Read the list of genomes
list_gff <- read_table(paste0(folder_data, "genomics/pangenome/old/list_gffs.txt"), col_names = F)
genome_ids <- list_gff$X1 %>% str_remove(".+/gff/") %>% str_remove(".gff")
length(genome_ids) # 31 genomes

# Read the list of core genes
folder_alignment <- paste0(folder_data, "genomics/pangenome/old/aligned_gene_sequences/")
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

write_csv(tibble(x = list_sccg), paste0(folder_data, "genomics_analysis/core_genes/list_sccg.csv"))
