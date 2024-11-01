#' This script computes the dxy using gene presence/absence as biallelic loci
#' Output one file per gradient
#' - per_acce_fst.csv which has the fst per-accessory-gene

library(tidyverse)
library(janitor)
#library(ape) # for reading multiple fasta
#library(poppr) # for reading snps into a genind object
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
make_longer_dist <- function (mdist) {
    mdist %>%
        as.matrix() %>%
        as_tibble() %>%
        mutate(genome_id1 = colnames(.)) %>%
        pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "dxy") %>%
        filter(genome_id1 < genome_id2)
}
compute_gcv_dxy <- function (tt) {
    #' This compute the dxy between the two populations
    temp <- tt$gpatl %>%
        left_join(isolates) %>%
        distinct(genome_id, .keep_all = T)

    ns <- table(temp$population)
    pop1 = names(ns)[1]
    pop2 = names(ns)[2]
    n1 = ns[1]
    n2 = ns[2]
    pop1_gid <- temp$genome_id[temp$population==pop1]
    pop2_gid <- temp$genome_id[temp$population==pop2]
    # number of gpa difference between sequences
    pi1 <- mean(dist(t(tt$gpa[,pop1_gid]), method = "manhattan")) # Within pop1 dxy
    pi2 <- mean(dist(t(tt$gpa[,pop2_gid]), method = "manhattan")) # Within pop2 dxy
    dall <- dist(t(tt$gpa[,-1]), method = "manhattan") # total
    dall_ln <- make_longer_dist(dall)
    # Between pops dxy
    dxy <- mean(as.matrix(dall)[1:n1,(n1+1):(n1+n2)])


    return(list(pop_dxy = tibble(pi1, pi2, dxy), ind_dxy = dall_ln))
}

set_name = "elev_med"
tt <- read_gpas(set_name)
dd <- compute_gcv_dxy(tt)
gcv_pop_dxy <- dd$pop_dxy
gcv_ind_dxy <- dd$ind_dxy

write_csv(gcv_pop_dxy, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gcv_pop_dxy.csv"))
write_csv(gcv_ind_dxy, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gcv_ind_dxy.csv"))


set_name <- "urbn_mel"
tt <- read_gpas(set_name)
dd <- compute_gcv_dxy(tt)
gcv_pop_dxy <- dd$pop_dxy
gcv_ind_dxy <- dd$ind_dxy

write_csv(gcv_pop_dxy, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gcv_pop_dxy.csv"))
write_csv(gcv_ind_dxy, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gcv_ind_dxy.csv"))



# for (set_name in c("elev_med", "urbn_mel")) {
#     #set_name = "elev_med"
#     #set_name = "urbn_mel"
#
#     tt <- read_gpas(set_name)
#     gcv_fst <- list()
#     for (gene in tt$gpa$gene) {
#         genind_object <- make_genind_from_gcv(tt$gpa, gene)
#         if (length(genind_object@loc.fac) == 1) {cat(gene, "is a core gene \n"); next}
#         ww <- wc(genind_object, diploid = F)
#         gcv_fst[[gene]] <- tibble(fst = ww$FST)
#         cat("\nProcessed:", gene)
#     }
#     per_acce_fst <- bind_rows(gcv_fst, .id = "gene") # per acce gene F_st results
#     write_csv(per_acce_fst, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/per_acce_fst.csv"))
#     # top_gene_or <- slice_top_genes(tidied_fisher)
#     # write_csv(top_gene_or, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/top_gene_or.csv"))
#
# }
#
#
