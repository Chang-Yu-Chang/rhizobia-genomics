#' This script computes the Fst using gene presence/absence as biallelic loci
#' Output one file per gradient
#' - per_acce_fst.csv which has the fst per-accessory-gene

library(tidyverse)
library(janitor)
library(ape) # for reading multiple fasta
library(poppr) # for reading snps into a genind object
library(hierfstat) # for computing Fst
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

make_genind_from_gcv <- function (gpa, gene_name) {
    gpa_i <- t(gpa[which(gpa$gene == gene_name),-1])
    pop <- isolates$population[match(rownames(gpa_i), isolates$genome_id)]
    df2genind(gpa_i, pop = pop, ploidy = 1)
}


for (set_name in c("elev_med", "urbn_mel")) {
    #set_name = "elev_med"
    #set_name = "urbn_mel"

    tt <- read_gpas(set_name)
    gcv_fst <- list()
    for (gene in tt$gpa$gene) {
        genind_object <- make_genind_from_gcv(tt$gpa, gene)
        if (length(genind_object@loc.fac) == 1) {cat(gene, "is a core gene \n"); next}
        ww <- wc(genind_object, diploid = F)
        gcv_fst[[gene]] <- tibble(fst = ww$FST)
        cat("\nProcessed:", gene)
    }
    per_acce_fst <- bind_rows(gcv_fst, .id = "gene") # per acce gene F_st results
    write_csv(per_acce_fst, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/per_acce_fst.csv"))
    # top_gene_or <- slice_top_genes(tidied_fisher)
    # write_csv(top_gene_or, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/top_gene_or.csv"))

}


