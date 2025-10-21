#'

library(tidyverse)
library(cowplot)
library(ape) # for reading alignment file
library(ggmsa) # for plotting msa
library(ggplotify) # for converting aplot into a grob
source(here::here("metadata.R"))

iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    arrange(contig_species)

list_genomes_symbiotic <- read_lines(paste0(folder_data, "genomics_analysis/count_snps/list_genomes_symbiotic.txt"))  # List of target genes
list_snps_symbiotic <- read_table(paste0(folder_data, "genomics_analysis/count_snps/list_snps_symbiotic.txt"))  # List of target genes
list_snps_symbiotic$GeneName <- str_remove(list_snps_symbiotic$GeneName, ".aln")
list_snps_symbiotic <- arrange(list_snps_symbiotic, GeneName)

n_nt_per_row = 500

for (gene_name in list_snps_symbiotic$GeneName) {
    # Read alignment
    aln <- read.FASTA(paste0(folder_genomics, "pangenome/aligned_gene_sequences/", gene_name, ".aln.fas"))
    names(aln) <- str_remove(names(aln), ";\\w+")
    aln <- aln[names(aln) %in% list_genomes_symbiotic] # subset those that are symbiotic
    # Re order by species
    genome_ids <- tibble(genome_id = names(aln)) %>%
        left_join(select(iso, genome_id, species = contig_species)) %>%
        arrange(species) %>%
        pull(genome_id)
    aln <- aln[genome_ids]

    #
    seq_length <- length(aln[[1]])
    n_rows <- ceiling(seq_length / n_nt_per_row)

    if (seq_length > 200) { # skip if the length of gene is shorter than 200
        p <- ggmsa(aln, start = 1, color = "Chemistry_NT", font = NULL, char_width = 0.5, seq_name = T) +
            facet_msa(field = n_nt_per_row) +
            theme(plot.background = element_rect(color = NA, fill = "white"))
        ggsave(paste0(folder_data, "genomics_analysis/count_snps/alignment_symbiotic/", gene_name, ".png"), p, width = 40, height = 4*n_rows)
        cat(gene_name)
    } else {
        cat(gene_name, " is shorter than ", n_nt_per_row)
        next
    }

}
