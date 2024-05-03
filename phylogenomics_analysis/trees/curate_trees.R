#' This script plots the consensus tree from iqtree

renv::load()
library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- isolates %>%
    left_join(isolates_contigs) %>%
    filter(!genome_id %in% c("g20", "g28"))

# Core gene tree
tr <- read.tree(paste0(folder_data, "genomics/mltree/isolates_core_b/aln.treefile"))
list_others <- c(paste0("g", c(20, 28, 38:43)), "em1022", "usda1106", "em1021", "wsm419")
tr <- tr %>% drop.tip(list_others)
tr <- root(tr, outgroup = "g15", resolve.root = TRUE)

# Compute gene counts
gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
gpatl <- gpat %>%
    pivot_longer(cols = -genome_id, names_to = "gene") %>%
    filter(value == 1) %>%
    filter(genome_id %in% tr$tip.label)
gene_order <- gpatl %>%
    group_by(gene) %>%
    dplyr::count() %>%
    arrange(desc(n)) %>%
    pull(gene)
gpatl <- gpatl %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    mutate(gene = factor(gene, gene_order))

# gpa tree
tr_gpa <- read.tree(paste0(folder_data, "genomics/mltree/isolates_gpa/aln.treefile"))
tr_gpa$tip.label <- gpat$genome_id[as.numeric(str_remove(tr_gpa$tip.label, "Seq"))]
tr_gpa <- root(tr_gpa, outgroup = "g15", resolve.root = TRUE)

save(tr, tr_gpa, gpatl, file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))






