#' This script cleans the gene presence/absence

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

# Gene presence-absence table
gpa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.Rtab")) %>%
    clean_names() # Rows are genes, columns are genomes
dim(gpa) # 26886 genes in union x 32-1 genomes
write_csv(gpa, paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))

# Transpose the gene presence-absence table
gpat <- gpa %>%
    pivot_longer(cols = -gene, names_to = "genome_id") %>%
    pivot_wider(names_from = gene, values_from = value, values_fill = 0)
dim(gpat) # 31 genomes x 26887 genes
write_csv(gpat, paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))

# Full directory for each gene
gpaf <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.csv")) %>%
    clean_names() %>%
    mutate(across(starts_with("g"), function (x) str_remove_all(x, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")))
dim(gpat) # 26886 genes x 34-3 genomes. 1) gene name 2) non unique gene name 3) annotation
gd <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_data.csv")) %>%
    clean_names() %>%
    mutate(annotation_id = str_remove_all(annotation_id, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")) %>%
    rename(genome_id = gff_file, contig_id )
dim(gd) # 225019 (gene x genome x contig) x 8 rows

# Structural variation
spa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/struct_presence_absence.Rtab")) %>%
    clean_names()
dim(spa) # 9775 strcutural variants x 32-1 genomes
write_csv(spa, paste0(folder_data, "genomics_analysis/gene_content/spa.csv"))

