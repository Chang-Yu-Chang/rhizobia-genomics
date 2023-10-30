#' This script

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. read data
list_g <- paste0("Chang_Q5C_", 1:19)
list_reads <- rep(list(NA), 19)

# 1. screen
i=1
mash_dist <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/screen.tab"), show_col_types = F, col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment"))
mash_dist %>%
    arrange(p_value)

list_mash <- rep(list(NA), length(list_g))


# 2. contigs
for (i in 1:length(list_g)) {
    i=1
    list_contigs <- list.files(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/"), pattern = "contig\\w+\\.tab") %>%
        str_replace(".tab", "")
    list_matches <- rep(list(NA), length(list_contigs))
    for (j in 1:length(list_contigs)) {
        list_matches[[j]] <- read_tsv(
            paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/", list_contigs[j],".tab"),
            show_col_types = F,
            col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment")) %>%
            arrange(p_value) %>%
            select(-median_multiplicity, -query_ID) %>%
            slice(1)
        #cat("\n", i)
    }
    list_mash[[i]] <- bind_rows(list_matches[nrow(list_matches) != 0], .id = "contig_id")
}

list_mash[-1]


#ggsave(paste0(folder_data, "temp/38-01-sabve.png"), p, width = 6, height = 5)
