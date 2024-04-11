#' This script assigns the taxonomy

renv::load()
library(tidyverse)
library(janitor)
library(seqinr)
library(ape)
source(here::here("metadata.R"))

# Read the list of genomes
list_gff <- read_table(paste0(folder_data, "genomics/pangenome/isolates/list_gffs.txt"), col_names = F)
genome_ids <- list_gff$X1 %>% str_remove(".+/gff/") %>% str_remove(".gff")
length(genome_ids) # 31 genomes

# Read the list of core genes
folder_alignment <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences/")
list_cg <- list.files(folder_alignment)

# Subset the core genes to only single-copy core genes
tb_core <- tibble(gene_name = list_cg, is_single_copy = NA)
for (i in 1:length(list_cg)) {
    aln <- read.alignment(paste0(folder_alignment, list_cg[i]), format = "fasta")
    list_seq_i <- str_extract(aln$nam, "\\w+;") %>% str_remove(";") %>% str_remove("_R_")
    tb_core$is_single_copy[i] <- length(unique(list_seq_i)) == length(list_seq_i)
    print(i)
}
sum(tb_core$is_single_copy) # 807 single-copy core genes
list_sccg <- list_cg[tb_core$is_single_copy]

write_csv(tibble(x = list_sccg), paste0(folder_data, "genomics_analysis/core_genes/list_sccg.csv"))

if (FALSE) {

# Generate a folder for single copy core genes
length(list_sccg)
source_dir <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences")
destination_dir <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences_single_copy")

if (!file.exists(destination_dir)) dir.create(destination_dir, recursive = TRUE)

# Loop through each file name in the vector
for (file_name in list_sccg) {
    # Construct the full path of the source file
    source_file <- file.path(source_dir, file_name)

    # Check if the file exists in the source directory
    if (file.exists(source_file)) {
        # Construct the full path of the destination file
        destination_file <- file.path(destination_dir, file_name)

        # Copy the file to the destination directory
        file.copy(source_file, destination_file, overwrite = TRUE)

        # Print a message indicating that the file has been copied
        cat("File", file_name, "has been copied to", destination_dir, "\n")
    } else {
        # Print a warning message if the file does not exist in the source directory
        cat("File", file_name, "does not exist in", source_dir, "\n")
    }
}

# Change the entry names

change_entry_names <- function(input_file, output_file) {
    # Read the input FASTA file
    fasta_content <- readLines(input_file)

    # Open the output file for writing
    output <- file(output_file, "w")

    # Loop through each line in the FASTA file
    for (line in fasta_content) {
        if (substr(line, 1, 1) == ">") {
            # If it's a header line, extract the part before the colon
            header_parts <- strsplit(line, ";")[[1]]
            new_header <- paste0("", header_parts[1], "\n")
            cat(new_header, file=output)
        } else {
            # Otherwise, write the line as is
            cat(line, "\n", file=output, append=TRUE)
        }
    }

    # Close the output file
    close(output)
}

# Example usage
# list_sccg[[1]]
# input_file <- paste0(destination_dir, "/", list_sccg[[1]])
# output_file <- paste0(destination_dir, "/", "test.aln.fas")
# change_entry_names(input_file, output_file)

for (alnfas in list_sccg) change_entry_names(paste0(destination_dir, "/", alnfas), paste0(destination_dir, "/", alnfas))









}
