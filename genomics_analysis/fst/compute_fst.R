#' This script computes Fst for a set of genes

library(tidyverse)
library(janitor)
library(ape) # for reading multiple fasta
library(poppr) # for reading snps into a genind object
library(mmod) # for computing Fst estimates: Nei's Gst, Hendrick's Gst, and Jost' D
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

read_aln <- function (aln_file) {
    paste0(aln_file) %>%
        read.dna(format = "fasta")
}
filter_snps <- function (genind_data, isolates) {
    #' this function removes snps with missing data
    # Clean index names
    indNames(genind_data) <- str_remove(indNames(genind_data), ";\\w+") %>% str_remove("_R_")
    # Assign population
    isolates_pop <- tibble(genome_id = indNames(genind_data)) %>% left_join(isolates, by = join_by(genome_id))
    pop(genind_data) <- isolates_pop$population
    # Remove SNP with missing data
    loci_to_keep <- apply(genind_data@tab, 2, function(x) !any(is.na(x)))
    genind_data_filtered <- genind_data[,loci_to_keep]

    return(genind_data_filtered)
}
fst_wrapper <- function (set_name) {
    #set_name <- "elev_med"
    #set_name <- "urbn_mel"
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    gene_list <- list_sccg$gene

    glen <- list()
    per_locus_fst_results <- list()
    per_gene_fst_results <- list()

    # Loop through each gene
    for (gene in gene_list) {
        #gene = "COQ5_4"
        #gene = "COQ3_2"
        # Read msa file
        alignment <- read_aln(paste0(folder_data, "genomics/pangenome/", set_name,"/aligned_gene_sequences/", gene, ".aln.fas"))
        glen[[gene]] <- ncol(alignment) # sequence length

        # Make a genind object
        gi <- DNAbin2genind(alignment)
        if (is.null(gi)) {cat(gene, "has no polymorphism. \n"); next}

        # Remove SNPs with missing data
        gi_filtered <- filter_snps(gi, isolates)
        if(ncol(gi_filtered@tab) == 0) {cat(gene, "has no polymorphism after removing samples with missing data. \n"); next}

        # Compute Fst
        #ww <- wc(gi_filtered, diploid = F) # Weir and Cockerham estimates of Fstatistics
        #hierfstat::boot.ppfst(gi_filtered)
        ww <- diff_stats(gi_filtered)

        per_locus_fst_results[[gene]] <- as_tibble(ww$per.locus) %>% mutate(location = rownames(ww$per.locus)) #extract_per_locus(ww, gi_filtered) # per locus F_st results
        per_gene_fst_results[[gene]] <- as_tibble(as.list(ww$global)) %>% mutate(n_snps = length(gi_filtered@loc.n.all)) # tibble(n_snps = length(gi_filtered@loc.n.all), fst = ww$FST) # gene-wide F_st results

        cat("\nProcessed:", gene, "\n")
    }

    gene_lengths <- unlist(glen) %>% enframe(name = "gene", value = "sequence_length")
    per_gene_fst <- bind_rows(lapply(names(per_gene_fst_results), function(gene) data.frame(gene = gene, per_gene_fst_results[[gene]]))) %>% as_tibble
    per_locus_fst <- bind_rows(lapply(names(per_locus_fst_results), function(gene) data.frame(gene = gene, per_locus_fst_results[[gene]]))) %>% as_tibble

    write_csv(gene_lengths, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_lengths.csv"))
    write_csv(per_gene_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/per_gene_fst.csv"))
    write_csv(per_locus_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/per_locus_fst.csv"))
}

fst_wrapper("elev_med")
fst_wrapper("urbn_mel")
