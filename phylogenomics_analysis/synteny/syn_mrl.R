#' This script follows Zhao 2021 pipeline for Syn-MRL. It includes four steps
#' 1. phylogenomic synteny network construction
#'  1a. all vs all reciprocal sequence similarity search for all annotated proteomes using DIAMOND
#'  1b. pairwise synteny block detection using MCScanX
#' 2. network clustering
#' 3. matrix representation
#' 4. ML based tree inference

renv::load()
library(tidyverse)
library(syntenet)
source(here::here("metadata.R"))

read_csd
gff <- read_delim(paste0(folder_data, "genomics/gff/genomes/g2.gff"), skip = 4, col_names = c("seq"))
gff <- read_gff3(paste0(folder_data, "genomics/gff/genomes/g2.gff"))
paste0(folder_data, "genomics/gff/genomes/g2.gff")
Biostrings::readAAStringSet()




if (F) {

gpa <- read_delim(paste0(folder_genomics, "pangenome/isolates/gene_presence_absence.Rtab"))
colnames(gpa)[-1]
m <- t(as.matrix(gpa[,-1]))
colnames(m) <- NULL
rownames(m) <- NULL

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
        sequence_name <- paste0("Seq", i)
        cat(sequence_name, " ", paste(binary_matrix[i, ], collapse = ""), "\n", file = con)
    }

    # Close the output file
    close(con)
}
dir.create(paste0(folder_genomics, "mltree/isolates_gpa"), showWarnings = F)
write_phylip(m, paste0(folder_genomics, "mltree/isolates_gpa/gpa_binary.phy"))
}

