#' This script

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. read data
list_g <- paste0("Chang_Q5C_", 1:19)
list_reads <- rep(list(NA), 19)

i=1
# Read screen of the sholw consensus
mash_dist <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/screen.tab"), show_col_types = F, col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment"))
mash_dist %>%
    arrange(p_value)

# Read contig_1
mash_dist <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/contig_2.tab"), show_col_types = F, col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment"))
mash_dist %>%
    #filter(p_value != 1) %>%
    arrange(p_value) %>%
    select(-median_multiplicity, -query_ID)

#mash_dist <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/contig_1.tab"), show_col_types = F, col_names = c("reference_id", "query_id", "mash_distance", "p_value", "matching_hashes"))

for (i in 1:19) list_reads[[i]] <- read_table(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "08-mash/distance.tab"), col_names = c("name", "phred", "length"))
raw_reads <- list_reads %>%
    bind_rows(.id = "genome_id") %>%
    mutate(genome_id = factor(genome_id, 1:20))

# subtract by 33
raw_reads$phred <- raw_reads$phred

raw_reads %>%
    filter(genome_id == 1) %>%
    pull(phred) %>%
    range()
