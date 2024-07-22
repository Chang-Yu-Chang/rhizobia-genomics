#' This script cleans the tree from iqtree output

renv::load()
library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- isolates %>%
    left_join(isolates_contigs) %>%
    filter(!genome_id %in% c("g20", "g28"))

# Core gene tree
tr <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/isolates_core_b/aln.treefile"))
list_others <- c(paste0("g", c(20, 28, 38:43)), "em1022", "usda1106", "em1021", "wsm419")
tr <- tr %>% drop.tip(list_others)
tr <- root(tr, outgroup = "g15", resolve.root = TRUE)

# gpa tree
tr_gpa <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/isolates_gpa/aln.treefile"))
tr_gpa$tip.label <- gpat$genome_id[as.numeric(str_remove(tr_gpa$tip.label, "Seq"))]
tr_gpa <- root(tr_gpa, outgroup = "g15", resolve.root = TRUE)

save(tr, tr_gpa, file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))






