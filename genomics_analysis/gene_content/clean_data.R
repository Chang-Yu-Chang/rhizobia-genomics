#' This script cleans the gene presence/absence and compute distance

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

gpa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.Rtab")) %>%
    clean_names() # Rows are genes, columns are genomes
write_csv(gpa, paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))

# Transpose the gene presence-absence table
gpat <- gpa %>%
    pivot_longer(cols = -gene, names_to = "genome_id") %>%
    pivot_wider(names_from = gene, values_from = value, values_fill = 0)

write_csv(gpat, paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))

# Full directory for each gene
gpaf <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.csv")) %>%
    clean_names() %>%
    mutate(across(starts_with("g"), function (x) str_remove_all(x, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")))
gd <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_data.csv")) %>%
    clean_names() %>%
    mutate(annotation_id = str_remove_all(annotation_id, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/"))

# gpaf %>%
#     select(gene, non_unique_gene_name, g2) %>%
#     drop_na(g2)
#
# gd %>%
#     filter(gff_file == "g2") %>%
#     filter(annotation_id == "g2.fasta_00001")

# Structural variation
spa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/struct_presence_absence.Rtab"))


