#' This script computes the dxy using gene presence/absence as biallelic loci
#' Output one file per gradient
#' - per_acce_fst.csv which has the fst per-accessory-gene

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
make_longer_dist <- function (mdist) {
    mdist %>%
        as.matrix() %>%
        as_tibble() %>%
        mutate(genome_id1 = colnames(.)) %>%
        pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "gcv_dxy") %>%
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
    gpaclw <- tt$gpacl %>%
        select(replicon_type, gene, genome_id) %>%
        mutate(temp = 1) %>%
        filter(replicon_type %in% c("chromosome", "pSymA", "pSymB", "pAcce")) %>%
        pivot_wider(names_from = genome_id, values_from = temp, values_fill = 0)

    # Whole genome
    ## number of gpa difference between sequences
    dall <- dist(t(tt$gpa[,-1]), method = "manhattan") # total
    dall_ln <- make_longer_dist(dall)
    ## Between pops dxy
    pi1 <- mean(dist(t(tt$gpa[,pop1_gid]), method = "manhattan")) # Within pop1 dxy
    pi2 <- mean(dist(t(tt$gpa[,pop2_gid]), method = "manhattan")) # Within pop2 dxy
    gcv_dxy <- mean(as.matrix(dall)[1:n1,(n1+1):(n1+n2)])

    # Replicon
    drep_ln <- gpaclw %>%
        nest(data = -replicon_type) %>%
        mutate(
            dall = map(data, ~dist(t(.x[,-1]), method = "manhattan")),
            dall_ln = map(dall, make_longer_dist)
        ) %>%
        unnest(dall_ln) %>%
        select(-data, -dall)

    # Number of accessory genes
    n_acce <- tt$gpacl %>%
        #distinct(replicon_type, gene) %>%
        group_by(replicon_type, gene) %>%
        count() %>%
        ungroup() %>%
        filter(n != max(n)) %>%
        group_by(replicon_type) %>%
        count(name = "n_accessory")

        #group_by(replicon_type) %>%

    return(list(pop_dxy = tibble(pi1, pi2, gcv_dxy), ind_dxy = dall_ln, rep_dxy = drep_ln, n_acce = n_acce))
}
gcv_dxy_wrapper <- function (set_name) {
    #set_name = "elev_med"
    #set_name <- "urbn_mel"
    tt <- read_gpas(set_name)
    dd <- compute_gcv_dxy(tt)

    write_csv(dd$pop_dxy, paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/pop_gcv_dxy.csv"))
    write_csv(dd$ind_dxy, paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/ind_gcv_dxy.csv"))
    write_csv(dd$rep_dxy, paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/rep_gcv_dxy.csv"))
    write_csv(dd$n_acce, paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/n_acce.csv"))

}

gcv_dxy_wrapper("elev_med")
gcv_dxy_wrapper("urbn_mel")

