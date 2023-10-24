#' This script

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. read data
list_g <- paste0(rep("Chang_Q5C_results", 19), "/Chang_Q5C_", 1:19, "/")
list_g[11] <- "Chang_Q5C_results_repeated/Chang_Q5C_11/"
list_g[18] <- "Chang_Q5C_results_repeated/Chang_Q5C_18/"
list_reads <- rep(list(NA), 19)

i=1
mash_dist <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "08-mash/distances.tab"), show_col_types = F, col_names = c("reference_id", "query_id", "mash_distance", "p_value", "matching_hashes"))
mash_dist

mash_dist %>%
    arrange(mash_distance, p_value) %>%
    slice(1:5)

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
