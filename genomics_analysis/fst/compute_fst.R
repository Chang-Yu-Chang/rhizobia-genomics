#' This script computes Fst for a set of genes

library(tidyverse)
library(apex) # for reading multiple fasta
library(poppr) # for processing fasta
library(mmod) # for computing Fst
library(pegas)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

set_name <- "urbn_mel"
list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
gene_list <- list_sccg$gene
per_locus_fst_results <- list()
gene_wide_fst_results <- list()

# Step 2: Loop through each gene
for (gene in gene_list) {
    # Construct file name - assuming the files are named as "gene.fasta"
    alignment_file <- paste0(folder_data, "genomics/pangenome/", set_name,"/aligned_gene_sequences/", gene, ".aln.fas")

    # Step 3: Check if the file exists
    if (file.exists(alignment_file)) {
        # Read the alignment file
        alignment <- read.dna(alignment_file, format = "fasta")

        # Convert to genind object
        genind_data <- DNAbin2genind(alignment)  # Adjust ploidy if necessary

        if (is.null(genind_data)) next; cat(gene, "has no polymorphism. \n")

        # Clean header names
        indNames(genind_data) <- str_remove(indNames(genind_data), ";\\w+") %>% str_remove("_R_")

        # Assign populations
        isolates_pop <- tibble(genome_id = indNames(genind_data)) %>% left_join(isolates, by = join_by(genome_id))
        pop(genind_data) <- isolates_pop$site_group

        # Compute per locus F_st
        fst <- diff_stats(genind_data)

        # Save per locus F_st results
        per_locus_fst_results[[gene]] <- fst$per.locus

        # Compute gene-wide F_st using overall allele frequencies
        gene_wide_fst <- fst$global

        # Save gene-wide F_st results
        gene_wide_fst_results[[gene]] <- gene_wide_fst

        # Optionally print progress
        cat("Processed:", gene, "\n")
    } else {
        cat("File not found for gene:", gene, "\n")
    }
}

per_locus_fst <- bind_rows(lapply(names(per_locus_fst_results), function(gene) data.frame(gene = gene, location = rownames(per_locus_fst_results[[gene]]), fst = per_locus_fst_results[[gene]]))) %>% as_tibble
gene_wide_fst <- bind_rows(lapply(names(gene_wide_fst_results), function(gene) data.frame(gene = gene, enframe(gene_wide_fst_results[[gene]])))) %>% as_tibble %>% pivot_wider

write_csv(per_locus_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/per_locus_fst.csv"))
write_csv(gene_wide_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_wide_fst.csv"))

