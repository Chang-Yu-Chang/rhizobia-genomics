#' This script binarizes the gene presence absence table purposed for iqtree input

renv::load()
library(tidyverse)
source(here::here("metadata.R"))

# Clean the gpaf names
gpaf <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpaf.csv"))
dim(gpaf) # 26886 x 34
gpafl <- gpaf %>%
    select(-non_unique_gene_name) %>%
    pivot_longer(-c(gene, annotation), names_to = "genome_id", values_to = "annotation_id", values_drop_na = T) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id) %>%
    mutate(annotation_id = str_remove(annotation_id, "_stop"))
dim(gpafl) # 213135

## Unpack genes with multiple copies in one genome
split_colons <- function(string) unlist(str_split(string, ";"))
gpafl_multi <- gpafl %>%
    filter(str_detect(annotation_id, ";")) %>%
    rowwise() %>%
    mutate(annotation_id = list(split_colons(annotation_id))) %>%
    unnest(annotation_id)
gpafl <- gpafl %>%
    filter(!str_detect(annotation_id, ";")) %>%
    bind_rows(gpafl_multi)
dim(gpafl) # 221,882



# Clean the gd names
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv")) %>% select(contig_id, replicon_type) %>% drop_na
gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gd.csv"))
dim(gd) # 225019 (gene x genome x contig) x 8 rows
gd <- gd %>% select(genome_id, contig_id, annotation_id)

# Join the two tables
gpacl <- gpafl %>%
    arrange(genome_id, annotation_id) %>%
    left_join(gd) %>%
    # Remove duplicate of multi-copy genes on one contig
    distinct(gene, contig_id) %>%
    # Filter for chromosome and two plasmids
    left_join(contigs)

length(unique(gpacl$genome_id))
length(unique(gpacl$contig_id))

# Matrix of contig genes
## chromosome
gpac_chrom <- gpacl %>%
    filter(replicon_type == "chromosome") %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = contig_id, values_from = value, values_fill = 0)

## pSymA
gpac_psyma <- gpacl %>%
    filter(replicon_type == "psymA like") %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = contig_id, values_from = value, values_fill = 0)

# pSymB
gpac_psymb <- gpacl %>%
    filter(replicon_type == "psymB like") %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = contig_id, values_from = value, values_fill = 0)



# Write a phylip files for each tree
write_phylip <- function(binary_matrix, output_file) {
    # Open the output file for writing
    con <- file(output_file, "w")

    # Get the number of sequences and sequence length
    num_sequences <- nrow(binary_matrix)
    seq_length <- ncol(binary_matrix)

    # Write the header
    cat(num_sequences, seq_length, "\n", file = con)

    # Write each sequence in PHYLIP format
    for (i in 1:num_sequences) {
        # PHYLIP sequence names are limited to 10 characters
        sequence_name <- rownames(m)[i]
        cat(sequence_name, " ", paste(binary_matrix[i, ], collapse = ""), "\n", file = con)
    }

    # Close the output file
    close(con)
}
wrapper_phylip <- function(bm, filename) {
    dir.create(paste0(folder_data, "phylogenomics_analysis/trees/mltree/", filename), showWarnings = F)
    write_phylip(bm, paste0(folder_data, "phylogenomics_analysis/trees/mltree/", filename, "/", filename, ".phy"))
}

# Genome gpa
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
m <- t(as.matrix(gpa[,-1]))
dim(m)
wrapper_phylip(m, "gpa_genomes")

# chromosome
m <- t(as.matrix(gpac_chrom[,c(-1,-2)]))
dim(m)
wrapper_phylip(m, "gpa_chrom")

# pSymA
m <- t(as.matrix(gpac_psyma[,c(-1,-2)]))
dim(m)
wrapper_phylip(m, "gpa_psyma")

# pSymB
m <- t(as.matrix(gpac_psymb[,c(-1,-2)]))
dim(m)
wrapper_phylip(m, "gpa_psymb")














