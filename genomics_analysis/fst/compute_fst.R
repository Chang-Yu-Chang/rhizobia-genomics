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

compute_fst <- function (af1, n1, af2, n2, aft) {
    #' Compute Fst from allele frequency
    ht = 1-(aft^2+(1-aft)^2) # total heterozygosity
    hs = n1/(n1+n2)*(1-af1^2-(1-af1)^2) + n2/(n1+n2)*(1-af2^2-(1-af2)^2) # within-population heterzygosity
    fst = 1-(hs/ht)
    return(tibble(ht = ht, hs = hs, fst = fst))

    # compute_fst(1/4, 4, 1/6, 6, 2/10)
    # ht = 1-(0.8^2+0.2^2)
    # hs = 0.4*(1-0.75^2-0.25^2) + 0.6*(1-(5/6)^2-(1/6)^2)
    # fst = 1-(hs/ht)
}

for (set_name in c("elev_med", "urbn_mel")) {
#set_name <- "elev_med"
#set_name <- "urbn_mel"
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
        suppressMessages(genind_data <- DNAbin2genind(alignment))

        if (is.null(genind_data)) {cat(gene, "has no polymorphism. \n"); next}

        # Clean header names
        indNames(genind_data) <- str_remove(indNames(genind_data), ";\\w+") %>% str_remove("_R_")

        # Assign populations
        isolates_pop <- tibble(genome_id = indNames(genind_data)) %>% left_join(isolates, by = join_by(genome_id))
        pop(genind_data) <- isolates_pop$site_group

        # Remove SNP with missing data
        loci_to_keep <- apply(genind_data@tab, 2, function(x) !any(is.na(x)))
        genind_data_filtered <- genind_data[,loci_to_keep]
        if(ncol(genind_data_filtered@tab) == 0) {cat(gene, "has no polymorphism after removing samples with missing data. \n"); next}

        # Compute per locus F_st
        snp_fst = list()
        for (snp in unique(genind_data_filtered@loc.fac)) {
            ii = str_detect(colnames(genind_data_filtered@tab), as.character(snp))
            ss <- genind_data_filtered@tab[,ii]
            pp <- tibble(population = pop(genind_data_filtered), allele1 = ss[,1], allele2 = ss[,2]) %>%
                mutate(allele1 = factor(allele1, c(0,1))) %>%
                group_by(population, allele1, .drop = F) %>%
                count() %>%
                group_by(population) %>%
                mutate(af = n/sum(n), n = sum(n)) %>%
                filter(allele1 == 1)

            snp_fst[[snp]] <- compute_fst(pp$af[1], pp$n[1], pp$af[2], pp$n[2], sum(ss[,1])/length(ss[,1]))
        }

        per_snp_fst <- bind_rows(snp_fst, .id = "location")
        # # Compute per locus F_st
        # fst <- diff_stats(genind_data_filtered)

        per_locus_fst_results[[gene]] <- per_snp_fst # per locus F_st results
        gene_wide_fst_results[[gene]] <- summarize(per_snp_fst, fst = mean(fst)) # gene-wide F_st results

        # Optionally print progress
        cat("\nProcessed:", gene, "\n")
    } else {
        cat("File not found for gene:", gene, "\n")
    }
}

#gene_wide_fst <- bind_rows(lapply(names(gene_wide_fst_results), function(gene) data.frame(gene = gene, enframe(gene_wide_fst_results[[gene]])))) %>% as_tibble %>% pivot_wider
gene_wide_fst <- bind_rows(lapply(names(gene_wide_fst_results), function(gene) data.frame(gene = gene, fst = gene_wide_fst_results[[gene]]))) %>% as_tibble
per_locus_fst <- bind_rows(lapply(names(per_locus_fst_results), function(gene) data.frame(gene = gene, per_locus_fst_results[[gene]]))) %>% as_tibble

write_csv(gene_wide_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_wide_fst.csv"))
write_csv(per_locus_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/per_locus_fst.csv"))



# Get the sequence length
glen <- list()
for (gene in gene_list) {
    alignment_file <- paste0(folder_data, "genomics/pangenome/", set_name,"/aligned_gene_sequences/", gene, ".aln.fas")
    sequences <- Biostrings::readDNAStringSet(alignment_file)
    glen[[gene]] <- length(sequences[[1]])
    cat("Processed:", gene, "\n")
}
gene_lengths <- unlist(glen) %>% enframe(name = "gene", value = "sequence_length")
write_csv(gene_lengths, paste0(folder_data, "genomics_analysis/fst/", set_name, "/gene_lengths.csv"))

}

#

