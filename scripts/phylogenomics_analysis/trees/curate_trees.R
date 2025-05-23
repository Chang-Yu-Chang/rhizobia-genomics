#' This script cleans the tree from iqtree output

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

curation_wrapper <- function () {
    tbtr <- tibble(tr_type = c("core", "gpa", "spa", "ani", "kmer"), tr = NA)
    tbtr$tr[1] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/combined_sccg/combined_sccg.treefile")))
    tbtr$tr[2] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/tr_gpa.tree")))
    tbtr$tr[3] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/tr_spa.tree")))
    tbtr$tr[4] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/tr_ani.tree")))
    tbtr$tr[5] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/tr_kmer.tree")))
    return(tbtr)
}

tbtr <- curation_wrapper()
save(tbtr, file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
