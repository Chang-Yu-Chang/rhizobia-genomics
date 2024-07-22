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
tr_seq_core <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/seq_core/seq_core.treefile"))
list_others <- c(paste0("g", c(20, 28, 38:43)), "em1022", "usda1106", "em1021", "wsm419")
tr_seq_core <- tr_seq_core %>% drop.tip(list_others)
tr_seq_core <- root(tr_seq_core, outgroup = "g15", resolve.root = TRUE)

# gpa trees
## genomes
tr_gpa_genomes <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_genomes/gpa_genomes.treefile"))
tr_gpa_genomes <- root(tr_gpa_genomes, outgroup = "g15", resolve.root = TRUE)

## chromosomes and plasmids
tr_gpa_chrom <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_chrom/gpa_chrom.treefile"))
tr_gpa_psyma <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_psyma/gpa_psyma.treefile"))
tr_gpa_psymb <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_psymb/gpa_psymb.treefile"))

#
save(tr_seq_core, tr_gpa_genomes, tr_gpa_chrom, tr_gpa_psyma, tr_gpa_psymb, file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))






if (F) {
    tr_gpa <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/isolates_gpa/aln.treefile"))
    tr_gpa$tip.label <- gpat$genome_id[as.numeric(str_remove(tr_gpa$tip.label, "Seq"))]
    tr_gpa <- root(tr_gpa, outgroup = "g15", resolve.root = TRUE)

}
