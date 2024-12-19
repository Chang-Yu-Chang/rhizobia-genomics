#' This script cleans the tree from iqtree output

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

curation_wrapper <- function (set_name) {
    tbtr <- tibble(set_name = rep(set_name, 4), tr_type = c("core", "gpa", "ani", "kmer"), tr = NA)
    tbtr$tr[1] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/combined_sccg/combined_sccg.treefile")))
    tbtr$tr[2] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_gpa.tree")))
    tbtr$tr[3] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_ani.tree")))
    tbtr$tr[4] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_kmer.tree")))
    return(tbtr)
}

tbtr1 <- curation_wrapper("elev_med")
tbtr2 <- curation_wrapper("urbn_mel")
tbtr <- bind_rows(tbtr1, tbtr2)
save(tbtr, file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
