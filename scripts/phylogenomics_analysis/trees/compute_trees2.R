#' Compute trees: gpa, ani, kmers

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
distl <- read_csv(paste0(folder_data, "genomics_analysis/distances/distl.csv"))


compute_gpa_tree <- function () {
    tt <- read_gpas()
    d <- dist(t(tt$gpa[,-1]))
    tr <- hclust(d, method = "ward.D2") %>% as.phylo()
    return(tr)
}
compute_spa_tree <- function () {
    tt <- read_gpas()
    d <- dist(t(tt$spa[,-1]))
    tr <- hclust(d, method = "ward.D2") %>% as.phylo()
    return(tr)
}
get_distl_gra <- function () {
    tt <- read_gpas()
    gids <- names(tt$gpa[,-1])
    distl %>%
        left_join(rename_all(isolates, ~paste0(.x, 1))) %>%
        left_join(rename_all(isolates, ~paste0(.x, 2))) %>%
        filter(genome_id1 %in% gids, genome_id2 %in% gids) %>%
        select(genome_id1, genome_id2, d_ani, d_kmer)
}
compute_ani_tree <- function (distl_gra) {
    d <- distl_gra %>%
        select(genome_id1, genome_id2, d_ani) %>%
        pivot_wider(names_from = genome_id2, values_from = d_ani)
    m <- as.dist(d[,-1])
    tr <- nj(m) %>% as.phylo()
    return(tr)
}
compute_kmer_tree <- function (distl_gra) {
    d <- distl_gra %>%
        select(genome_id1, genome_id2, d_kmer) %>%
        pivot_wider(names_from = genome_id2, values_from = d_kmer)
    m <- as.dist(d[,-1])
    tr <- nj(m) %>% as.phylo()
    return(tr)
}


tr_gpa <- compute_gpa_tree()
tr_spa <- compute_spa_tree()

distl_gra <- get_distl_gra()
tr_ani <- compute_ani_tree(distl_gra)
tr_kmer <- compute_kmer_tree(distl_gra)
write.tree(tr_gpa, paste0(folder_data, "phylogenomics_analysis/trees/tr_gpa.tree"))
write.tree(tr_spa, paste0(folder_data, "phylogenomics_analysis/trees/tr_spa.tree"))
write.tree(tr_ani, paste0(folder_data, "phylogenomics_analysis/trees/tr_ani.tree"))
write.tree(tr_kmer, paste0(folder_data, "phylogenomics_analysis/trees/tr_kmer.tree"))
