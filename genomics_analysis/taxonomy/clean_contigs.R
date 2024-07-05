#' This script joins contig data

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

# Replicon size by nt
genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))

# Contig level taxonomy
b_genome <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/b_genome.csv")) %>%
    rename(contig_id = qseqid) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id)) %>%
    select(genome_id, contig_id, species, strain, replicon, pident, length)

# Gene content
gpaf <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.csv")) %>%
    clean_names() %>%
    mutate(across(starts_with("g"), function (x) str_remove_all(x, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")))
dim(gpaf) # 26886 genes in the pangenome
gd <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_data.csv")) %>%
    clean_names() %>%
    mutate(annotation_id = str_remove_all(annotation_id, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")) %>%
    rename(genome_id = gff_file, contig_id = scaffold_name) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id))
dim(gd) # 225019 genes across genomes and which contig they are from. The longer format of gpaf

# Number of genes
gdn <- gd %>%
    group_by(genome_id, contig_id) %>%
    count(name = "n_genes")

# Join the datasets
contigs <- genomes %>%
    left_join(gdn) %>%
    left_join(b_genome) %>%
    drop_na(species, n_genes) %>%
    mutate(frac_blas = round(length / contig_length * 100, 2))

write_csv(contigs, paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))

# Join the by contig genes
gpafl <- gpaf %>%
    mutate(gene_id = 1:n()) %>%
    select(gene_id, gene, starts_with("g")) %>%
    pivot_longer(cols = -c(gene_id, gene), names_to = "genome_id", values_to = "annotation_id", values_drop_na = T) %>%
    separate_rows(annotation_id, sep = ";")

gpafl %>%
    filter(gene_id==4)

gpaflc <- gpafl %>%
    left_join(select(gd, genome_id, contig_id, annotation_id))
    # Somehow some gene are not matched?
    drop_na(contig_id)

gpaflc %>%
    select(-gene) %>%
    group_by(contig_id) %>%
    count() %>%
    view
    #select(gene_id, gene, contig_id, annotation_id) %>%
    #mutate(temp = 1) %>%
    pivot_wider(names_from = contig_id, values_from =  annotation_id)
















