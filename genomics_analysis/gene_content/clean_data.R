#' This script cleans the gene presence/absence from panaroo output

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

# Compute gene counts
gpatl <- gpat %>%
    pivot_longer(cols = -genome_id, names_to = "gene") %>%
    filter(value == 1)
gene_order <- gpatl %>%
    group_by(gene) %>%
    dplyr::count() %>%
    arrange(desc(n)) %>%
    pull(gene)
gpatl <- gpatl %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    mutate(gene = factor(gene, gene_order))
dim(gpatl) # 213135 3
length(gene_order) # 26886 genes order by copy number
write_csv(gpatl, paste0(folder_data, "genomics_analysis/gene_content/gpatl.csv"))
write_csv(tibble(gene = gene_order), paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv"))

# Full directory for each gene
gpaf <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.csv")) %>%
    clean_names() %>%
    mutate(across(starts_with("g"), function (x) str_remove_all(x, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")))
dim(gpaf) # 26886 genes x 34 genomes. 1) gene name 2) non unique gene name 3) annotation
write_csv(gpaf, paste0(folder_data, "genomics_analysis/gene_content/gpaf.csv"))

# Genes and the contigs they are from
gd <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_data.csv")) %>%
    clean_names() %>%
    mutate(annotation_id = str_remove_all(annotation_id, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")) %>%
    rename(genome_id = gff_file, contig_id = scaffold_name)
dim(gd) # 225019 (gene x genome x contig) x 8 rows
write_csv(gd, paste0(folder_data, "genomics_analysis/gene_content/gd.csv"))

# Structural variation
spa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/struct_presence_absence.Rtab")) %>%
    clean_names()
dim(spa) # 9775 structural variants x 32-1 genomes
write_csv(spa, paste0(folder_data, "genomics_analysis/gene_content/spa.csv"))











