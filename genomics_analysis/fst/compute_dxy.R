#' This script computes Fst for a set of genes

library(tidyverse)
library(StAMPP)
library(apex) # for reading multiple fasta
library(ape)  # for genetic distance
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


# Example hypothetical data
# Create a genotype matrix (rows are SNPs, columns are individuals)
# 0 = homozygous for allele 1, 1 = heterozygous, 2 = homozygous for allele 2
genotypes <- matrix(c(0, 1, 1, 0, 2, 0, 2, 1, 0, 1, 1, 0), nrow = 6, ncol = 2)

# Assign populations (First 3 for population A, last 3 for population B)
pop_labels <- c(rep("A", 3), rep("B", 3))

# Create a genind object, if using adegenet package:
library(adegenet)
#genind_obj <- df2genind(genotypes, ploidy=2, ind.names=1: ncol(genotypes), loc.names=1:nrow(genotypes))

# For this example, calculate Dxy manually
calculate_dxy <- function(geno_matrix, popA, popB) {
    n_sites <- nrow(geno_matrix)
    n_A <- length(popA)
    n_B <- length(popB)

    total_dxy <- 0

    for (i in 1:n_sites) {
        # Get the alleles for this SNP in each population
        alleles_A <- geno_matrix[i, popA]
        alleles_B <- geno_matrix[i, popB]

        # Calculate the number of differences
        dxy <- sum(outer(alleles_A, alleles_B, "!=")) / (n_A * n_B)
        total_dxy <- total_dxy + dxy
    }

    # Average Dxy across sites
    normalized_dxy <- total_dxy / n_sites
    return(normalized_dxy)
}

# Population indices (assuming columns 1-3 belong to Pop A and columns 4-6 to Pop B)
popA_indices <- 1:3
popB_indices <- 4:6

# Calculate Dxy
dxy_result <- calculate_dxy(genotypes, popA_indices, popB_indices)
print(paste("Dxy:", dxy_result))



for (set_name in c("elev_med", "urbn_mel")) {
    set_name <- "elev_med"
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

            # Compute per locus nucleotide difference
            snp_fst = list()
            for (snp in unique(genind_data_filtered@loc.fac)) {
                ii = str_detect(colnames(genind_data_filtered@tab), as.character(snp))
                ss <- genind_data_filtered@tab[,ii]
                pp <- tibble(population = pop(genind_data_filtered), allele1 = ss[,1], allele2 = ss[,2]) %>%
                    mutate(allele1 = factor(allele1, c(0,1)))



                snp_fst[[snp]] <- compute_dxy(pp$af[1], pp$n[1], pp$af[2], pp$n[2], sum(ss[,1])/length(ss[,1]))
            }
        }
    }
}



# import genotype data and convert to allele frequecies
data(potato.mini, package="StAMPP")
potato.mini <- as_tibble(potato.mini)[,1:10]
potato.freq <- stamppConvert(potato.mini, "r")
dim(potato.freq)
class(potato.mini)
# Calculate genetic distance between individuals
potato.D.ind <- stamppNeisD(potato.freq, FALSE)
# Calculate genetic distance between populations
potato.D.pop <- stamppNeisD(potato.freq, TRUE)
