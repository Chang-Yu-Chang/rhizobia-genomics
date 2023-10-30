#' Plot the sourmash results

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))


# 0. Read data
list_sm <- rep(list(NA), 19) # list of sourmash
for (i in 2:length(list_sm)) {
    list_sm[[i]] <- read_csv(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_", i, "/09-sourmash/gathered.csv"), show_col_types = F) %>%
        mutate(query_md5 = as.character(query_md5)) %>%
        mutate(genome_id = i)
}

tb_sm <- bind_rows(list_sm[-1]) %>%  # tibble of sourmash
    select(genome_id, everything())


# 1. summary
tb_sm %>%
    group_by(genome_id) %>%
    summarize(n_matches = n())


#
tb_sm_tax <- tb_sm %>%
    group_by(genome_id) %>%
    slice(1) %>%
    #select(genome_id, intersect_bp, name) %>%
    mutate(`intersect_bp (Mbp)` = intersect_bp/1000000) %>%
    mutate(name = str_replace(name, "GC[A|F]_\\d+\\.1 ", "")) %>%
    #mutate(contig_id = 1:n()) %>%
    select(genome_id, `intersect_bp (Mbp)`, f_orig_query, f_match, name)

write_csv(tb_sm_tax, paste0(folder_data, "temp/39-tb_sm_tax.csv"))
