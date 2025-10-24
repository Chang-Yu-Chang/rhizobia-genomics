#' Assign the gpa file to replicon based on which contig its on

library(tidyverse)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
contigs <- read_csv(paste0(folder_genomics, "ibd/contigs.csv"))
gpar <- read_csv(file.path(folder_genomics, "pangenome/gene_content/gpar.csv"), show_col_types = FALSE)
gd <- read_csv(file.path(folder_genomics, "pangenome/gene_content/gd.csv"), show_col_types = FALSE)

# prepare
gd <- mutate(gd, annotation_id = basename(annotation_id))
contigs <- mutate(contigs, contig_id = paste0(genome_id, "_", qseqid))

gpafl <- gpar %>%
    select(gene, matches("g\\d+")) %>%
    pivot_longer(cols = matches("g\\d+"), names_to = "genome_id", values_to = "annotation_id", values_drop_na = TRUE) %>%
    # split cells with multiple entries separated by ";"
    separate_rows(annotation_id, sep = ";") %>%
    mutate(annotation_id = str_remove(annotation_id, "_stop"), annotation_id = basename(annotation_id))

#
gpacl <- gpafl %>%
    left_join(gd %>% select(genome_id, contig_id, annotation_id),  by = c("genome_id", "annotation_id")) %>%
    distinct(gene, contig_id, .keep_all = TRUE) %>%
    left_join(contigs, by = c("genome_id", "contig_id"))

gpacl <- write_csv(gpacl, paste0(folder_genomics, "ibd/gpacl.csv"))
