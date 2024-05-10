#' This script binarizes the gene presence absence table purposed for iqtree input

renv::load()
library(tidyverse)
source(here::here("metadata.R"))

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
