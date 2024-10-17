#' This script computes Fst for a set of genes

library(tidyverse)
library(poppr) # for reading snps into a genind object
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
compute_af <- function (gi_filtered, snp) {
    #' This function computes the allele frequency per population
    gi_temp <- gi_filtered[, str_detect(colnames(gi_filtered@tab), snp)]
    ns <- table(gi_temp@pop)
    pop1 = names(ns)[1]
    pop2 = names(ns)[2]
    n1 = ns[1]
    n2 = ns[2]
    a1 = gi_temp@tab[which(gi_temp@pop == pop1),1]
    af1 = mean(a1)
    a2 = gi_temp@tab[which(gi_temp@pop == pop2),1]
    af2 = mean(a2)
    aft = mean(gi_temp@tab[,1])

    return(tibble(af1, n1, af2, n2, aft))
}
compute_fst <- function (af1, n1, af2, n2, aft, fpc = T) {
    #' Compute Fst from allele frequency. Same as poppr(gi)$Hexp

    # Correct for finite/small population size
    if (fpc) {
        h1 = n1/(n1-1)*(1-af1^2-(1-af1)^2) # pop 1
        h2 = n2/(n2-1)*(1-af2^2-(1-af2)^2) # pop 2
        ht = (n1+n2)/(n1+n2-1)*(1-aft^2-(1-aft)^2) # total heterozygosity
    } else if (fpc == F){
        h1 = 1-af1^2-(1-af1)^2 # pop 1
        h2 = 1-af2^2-(1-af2)^2 # pop 2
        ht = 1-aft^2-(1-aft)^2 # total heterozygosity
    }
    hs = (n1*h1+n2*h2)/(n1+n2) # within-population heterozygosity
    fst = 1-(hs/ht)
    return(tibble(h1, h2, ht, hs, fst))
}

for (set_name in c("elev_med", "urbn_mel")) {
    #set_name <- "elev_med"
    #set_name <- "urbn_mel"
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    gene_list <- list_sccg$gene

    glen <- list()
    per_locus_fst_results <- list()
    gene_wide_fst_results <- list()

    # Loop through each gene
    for (gene in gene_list) {

        # Read msa file
        alignment <- read_aln(paste0(folder_data, "genomics/pangenome/", set_name,"/aligned_gene_sequences/", gene, ".aln.fas"))
        glen[[gene]] <- ncol(alignment) # sequence length

        # Make a genind object
        gi <- DNAbin2genind(alignment)
        if (is.null(gi)) {cat(gene, "has no polymorphism. \n"); next}

        # Remove SNPs with missing data
        gi_filtered <- filter_snps(gi, isolates)
        if(ncol(gi_filtered@tab) == 0) {cat(gene, "has no polymorphism after removing samples with missing data. \n"); next}

        # Compute per locus Fst
        list_snps <- names(gi_filtered@all.names)
        snp_fst <- list()
        for (snp in list_snps) {
            afs <- compute_af(gi_filtered, snp)
            snp_fst[[snp]] <- compute_fst(afs$af1, afs$n1, afs$af2, afs$n2, afs$aft)
        }
        per_snp_fst <- bind_rows(snp_fst, .id = "location")

        per_locus_fst_results[[gene]] <- per_snp_fst # per locus F_st results
        gene_wide_fst_results[[gene]] <- summarize(per_snp_fst, n_snps = n(), fst = mean(fst)) # gene-wide F_st results

        cat("\nProcessed:", gene, "\n")
    }

    gene_lengths <- unlist(glen) %>% enframe(name = "gene", value = "sequence_length")
    gene_wide_fst <- bind_rows(lapply(names(gene_wide_fst_results), function(gene) data.frame(gene = gene, fst = gene_wide_fst_results[[gene]]))) %>% as_tibble
    per_locus_fst <- bind_rows(lapply(names(per_locus_fst_results), function(gene) data.frame(gene = gene, per_locus_fst_results[[gene]]))) %>% as_tibble

    write_csv(gene_lengths, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_lengths.csv"))
    write_csv(gene_wide_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_wide_fst.csv"))
    write_csv(per_locus_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/per_locus_fst.csv"))
}
