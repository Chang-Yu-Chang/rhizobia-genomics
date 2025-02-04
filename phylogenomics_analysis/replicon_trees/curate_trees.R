#' This script cleans the tree from iqtree output

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

curation_wrapper <- function (set_name) {
    tbtr <- tibble(set_name = rep(set_name, 8), tr_type = c("core", "gpa"), replicon = c("chromosome", "pSymA", "pSymB", "pAcce"), tr = NA)
    tbtr$tr[1] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/combined_sccg/combined_sccg.treefile")))
    tbtr$tr[2] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_gpa.tree")))
    return(tbtr)
}

set_name = "elev_med"
tbtr1 <- tibble(set_name = rep(set_name, 8), tr_type = rep(c("core", "gpa"), each = 4), replicon_type = rep(c("chromosome", "pSymA", "pSymB", "pAcce"), 2), tr = NA)
for (i in 1:4) tbtr1$tr[i] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/", tbtr1$replicon_type[i], "/combined_sccg/combined_sccg.treefile")))
for (i in 5:8) tbtr1$tr[i] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/", tbtr1$replicon_type[i], "/tr_gpa.tree")))


set_name = "urbn_mel"
tbtr2 <- tibble(set_name = rep(set_name, 8), tr_type = rep(c("core", "gpa"), each = 4), replicon_type = rep(c("chromosome", "pSymA", "pSymB", "pAcce"), 2), tr = NA)
for (i in 1:3) tbtr2$tr[i] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/", tbtr2$replicon_type[i], "/combined_sccg/combined_sccg.treefile")))
for (i in 5:7) tbtr2$tr[i] <- list(read.tree(paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/", tbtr2$replicon_type[i], "/tr_gpa.tree")))
tbtr2 <- filter(tbtr2, !(replicon_type == "pAcce"))

tbtr <- bind_rows(tbtr1, tbtr2)
save(tbtr, file = paste0(folder_data, "phylogenomics_analysis/replicon_trees/trees.rdata"))
