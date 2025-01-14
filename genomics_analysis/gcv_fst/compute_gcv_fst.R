#' This script computes Fst for a set of genes

library(tidyverse)
library(janitor)
library(ape) # for reading multiple fasta
library(poppr) # for reading snps into a genind object
library(mmod) # for computing Fst estimates: Nei's Gst, Hendrick's Gst, and Jost' D
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))


gcv_fst_wrapper <- function (set_name) {
    #set_name <- "elev_med"
    #set_name <- "urbn_mel"

    tt <- read_gpas(set_name)

    # Make genind object from table
    gi <- df2genind(t(tt$gpa[,-1]), ploidy = 1)
    isolates_pop <- tibble(genome_id = indNames(gi)) %>% left_join(isolates, by = join_by(genome_id))
    pop(gi) <- isolates_pop$population
    ww <- diff_stats(gi)
    per_acce_fst <- as_tibble(ww$per.locus) %>% mutate(loc_id = rownames(ww$per.locus), gene = tt$gpa$gene)
    per_genome_fst <- as_tibble(as.list(ww$global))
    write_csv(per_acce_fst, paste0(folder_data, "genomics_analysis/gcv_fst/", set_name, "/per_acce_fst.csv"))
    write_csv(per_genome_fst, paste0(folder_data, "genomics_analysis/gcv_fst/", set_name, "/per_genome_fst.csv"))

    #gene_lengths <- unlist(glen) %>% enframe(name = "gene", value = "sequence_length")
    #per_gene_fst <- bind_rows(lapply(names(per_gene_fst_results), function(gene) data.frame(gene = gene, per_gene_fst_results[[gene]]))) %>% as_tibble
    # per_locus_fst <- bind_rows(lapply(names(per_locus_fst_results), function(gene) data.frame(gene = gene, per_locus_fst_results[[gene]]))) %>% as_tibble
    # write_csv(gene_lengths, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_lengths.csv"))
    # write_csv(per_gene_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/per_gene_fst.csv"))
    # write_csv(per_locus_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/per_locus_fst.csv"))
}

gcv_fst_wrapper("elev_med")
gcv_fst_wrapper("urbn_mel")
