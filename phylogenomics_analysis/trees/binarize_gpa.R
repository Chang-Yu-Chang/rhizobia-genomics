#' This script binarizes the gene presence absence table purposed for iqtree input

renv::load()
library(tidyverse)
source(here::here("metadata.R"))

contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv")) %>% select(contig_id, replicon_type) %>% drop_na
gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpacl.csv")) %>%
    filter(contig_id %in% contigs$contig_id)
dim(gpacl)


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

# Structural variants
spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/spa.csv"))
m <- t(as.matrix(spa[,-1]))
dim(m)
wrapper_phylip(m, "spa_genomes")














