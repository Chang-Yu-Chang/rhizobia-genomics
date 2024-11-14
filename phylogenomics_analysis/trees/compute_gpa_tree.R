#' Compute GCV tree

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- isolates %>% left_join(isolates_contigs)


compute_gpa_tree <- function (set_name) {
    tt <- read_gpas(set_name)
    d <- dist(t(tt$gpa[,-1]))
    tr <- hclust(d, method = "ward.D2") %>% as.phylo()
    return(tr)

}

set_name <- "elev_med"
tr <- compute_gpa_tree(set_name)
write.tree(tr, paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_gpa.tree"))

set_name <- "urbn_mel"
tr <- compute_gpa_tree(set_name)
write.tree(tr, paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_gpa.tree"))

