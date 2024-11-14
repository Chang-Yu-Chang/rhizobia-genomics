#' This script cleans the tree from iqtree output

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- isolates %>% left_join(isolates_contigs)


set_name <- "elev_med"
tr_elev_med_core <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/combined_sccg/combined_sccg.treefile"))
tr_elev_med_gpa <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_gpa.tree"))

set_name <- "urbn_mel"
tr_urbn_mel_core <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/combined_sccg/combined_sccg.treefile"))
tr_urbn_mel_gpa <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/", set_name, "/tr_gpa.tree"))

save(
    tr_elev_med_core,
    tr_elev_med_gpa,
    tr_urbn_mel_core,
    tr_urbn_mel_gpa,
    file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata")
)

if (F) {

# Core gene tree

# tr_seq_core <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/seq_core/seq_core.treefile"))
#tr_seq_core <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/deprecated/isolates_core_b/aln.treefile"))
# list_others <- c(paste0("g", c(20, 28, 38:43)), "em1022", "usda1106", "em1021", "wsm419")
# tr_seq_core <- tr_seq_core %>% drop.tip(list_others)
# tr_seq_core <- root(tr_seq_core, outgroup = "g15", resolve.root = TRUE)

# gpa trees
## genomes
tr_gpa_genomes <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_genomes/gpa_genomes.treefile"))
tr_gpa_genomes <- root(tr_gpa_genomes, outgroup = "g15", resolve.root = TRUE)

## chromosomes and plasmids
tr_gpa_chrom <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_chrom/gpa_chrom.treefile"))
tr_gpa_psyma <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_psyma/gpa_psyma.treefile"))
tr_gpa_psymb <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/gpa_psymb/gpa_psymb.treefile"))

# structural variants presence and absence
tr_spa_genomes <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/spa_genomes/spa_genomes.treefile"))

#
save(tr_seq_core, tr_gpa_genomes, tr_gpa_chrom, tr_gpa_psyma, tr_gpa_psymb, tr_spa_genomes,
     file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

}
