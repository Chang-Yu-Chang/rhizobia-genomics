#' This script computes Fst for a set of genes

library(tidyverse)
library(poppr) # for reading snps into a genind object and computing dxy
library(ape) # for reading multiple fasta
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))

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
make_longer_dist <- function (mdist) {
    mdist %>%
        as.matrix() %>%
        as_tibble() %>%
        mutate(genome_id1 = colnames(.)) %>%
        pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "dxy") %>%
        filter(genome_id1 < genome_id2)
}
compute_dxy <- function (gi_filtered) {
    #' This compute the dxy between the two populations
    ns <- table(gi_filtered@pop)
    pop1 = names(ns)[1]
    pop2 = names(ns)[2]
    n1 = ns[1]
    n2 = ns[2]

    # number of nucleotide difference between sequences
    pi1 <- mean(diss.dist(popsub(gi_filtered, pop1))) # Within pop1 dxy
    pi2 <- mean(diss.dist(popsub(gi_filtered, pop2))) # Within pop2 dxy
    dall <- bitwise.dist(gi_filtered, percent = F)
    # mean(dall) should be the same as the sum of total heterozygosity across site sum(per_snp_fst$ht)
    # Individual sample dxy
    dall_ln <- make_longer_dist(dall)
    # Between pops dxy
    dxy <- mean(as.matrix(dall)[1:n1,(n1+1):(n1+n2)])

    return(list(pop_dxy = tibble(pi1, pi2, dxy), ind_dxy = dall_ln))
}
ibd_wrapper <- function (set_name) {
    #set_name <- "elev_med"
    #set_name <- "urbn_mel"
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    gene_list <- list_sccg$gene
    pop_dxy_results <- list()
    ind_dxy_results <- list()

    for (gene in gene_list) {
        # Read msa file
        alignment <- read_aln(paste0(folder_data, "genomics/pangenome/", set_name,"/aligned_gene_sequences/", gene, ".aln.fas"))

        # Make a genind object
        gi <- DNAbin2genind(alignment)
        if (is.null(gi)) {cat(gene, "has no polymorphism. \n"); next}

        # Remove SNPs with missing data
        gi_filtered <- filter_snps(gi, isolates)
        if(ncol(gi_filtered@tab) == 0) {cat(gene, "has no polymorphism after removing samples with missing data. \n"); next}

        # Compute pi and dxy
        dd <- compute_dxy(gi_filtered)
        pop_dxy_results[[gene]] <- dd$pop_dxy
        ind_dxy_results[[gene]] <- dd$ind_dxy

        cat("\nProcessed:", gene, "\n")
    }

    gene_pop_dxy <- bind_rows(pop_dxy_results, .id = "gene")
    gene_ind_dxy <- bind_rows(ind_dxy_results, .id = "gene")
    write_csv(gene_pop_dxy, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_pop_dxy.csv"))
    write_csv(gene_ind_dxy, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_ind_dxy.csv"))

}

ibd_wrapper("elev_med")
ibd_wrapper("urbn_mel")

