#' This script cleans the tree from iqtree output

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

tbtr <- tibble(tr_type = c("core", "gpa", "spa"), tr = NA)
tbtr$tr[1] <- list(read.tree(paste0(folder_genomics, "pangenome/trees/combined_sccg/combined_sccg.treefile")))
tbtr$tr[2] <- list(read.tree(paste0(folder_genomics, "pangenome/trees/tr_gpa.tree")))
tbtr$tr[3] <- list(read.tree(paste0(folder_genomics, "pangenome/trees/tr_spa.tree")))

save(tbtr, file = paste0(folder_genomics, "pangenome/trees/trees.rdata"))
