#' This script

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 0. read data
isolates <- read_csv(paste0(folder_data, "temp/42-isolates.csv"), show_col_types = F) %>%
    mutate(genome_name = str_replace(genome_id, "g", "Chang_Q5C_") %>% factor(paste0("Chang_Q5C_", 1:20)))

# 0.1 read the whole genome mash (taking all contigs of a genome)
#' Aggregate the mash screen result
#' This only needs to run once
list_g <- paste0("Chang_Q5C_", 1:19)
list_screen <- rep(list(NA), 19)

for (i in 1:19) {
    mash_dist <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/screen.tab"), show_col_types = F, col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment"))
    list_screen[[i]] <- mash_dist %>%
        arrange(p_value) %>%
        mutate(genome_id = paste0("g", i)) %>%
        select(genome_id, everything())
}
mash_g <- list_screen %>% bind_rows()
write_csv(mash_g, paste0(folder_data, "temp/38-mash_g.csv"))

mash_g <- read_csv(paste0(folder_data, "temp/38-mash_g.csv"), show_col_types = F)
isolates_rhizo <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F)

# Clean up
mash_g <- mash_g %>%
    left_join(isolates_rhizo) %>%
    filter(genus == "Ensifer") %>%
    mutate(genome_id = factor(genome_id, paste0("g", 1:19)))

# 0.2 read the mash outcome for individual contigs
#' Aggregate the mash screnn result
#' This only needs to run once
list_g <- paste0("Chang_Q5C_", 1:19)
list_screen <- rep(list(NA), 19)

for (i in 1:19) {
    list_contig_names <- list.files(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/08-mash/"), pattern = "contig", full.names = T) %>%
        str_subset(".tab") %>% unique()
    list_contig_screen <- rep(list(NA), length(list_contig_names))
    for (j in 1:length(list_contig_names)) {
        mash_dist <- read_tsv(list_contig_names[j], show_col_types = F, col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment"))
        if (nrow(mash_dist)==0) next
        list_contig_screen[[j]] <- mash_dist %>%
            arrange(p_value) %>%
            mutate(contig_id = paste0("contig_", j)) %>%
            select(contig_id, everything())
    }
    list_contig_screen <- list_contig_screen[!is.na(list_contig_screen)]
    list_screen[[i]] <- bind_rows(list_contig_screen) %>%
        mutate(genome_id = paste0("g", i)) %>%
        select(genome_id, contig_id, everything())
}
mash_c <- list_screen %>% bind_rows()
write_csv(mash_c, paste0(folder_data, "temp/38-mash_c.csv"))

# 1. genome top hits
mash_g_top <- mash_g %>%
    group_by(genome_id) %>%
    slice(1:3)
write_csv(mash_g_top, paste0(folder_data, "temp/38-mash_g_top.csv"))


# 2. contig top hits
mash_c_top <- mash_c %>%
    group_by(genome_id, contig_id) %>%
    slice(1:3)
write_csv(mash_c_top, paste0(folder_data, "temp/38-mash_c_top.csv"))


# 3. find the contigs with top length
mash_c_top <- read_csv(paste0(folder_data, "temp/38-mash_c_top.csv"), show_col_types = F)

egcct <- read_csv(paste0(folder_data, "temp/45-egcct.csv"), show_col_types = F)
egcct_c <- egcct %>%
    distinct(genome_id, contig, contig_ordered, contig_type)

mash_c_top %>%
    select(genome_id, contig = contig_id, identity, query_comment) %>%
    left_join(egcct_c) %>%
    mutate(genome_id = factor(genome_id, paste0("g", 1:20))) %>%
    arrange(genome_id, contig_type) %>%
    drop_na() %>%
    #filter(str_detect(query_comment, "Ensifer|Sinorhizobium")) %>%
    select(genome_id, contig_type, query_comment) %>%
    group_by(genome_id, contig_type) %>%
    mutate(mash_hits = paste0("mash", 1:n())) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(genome_id, contig_type), names_from = mash_hits, values_from = query_comment) %>%
    left_join(select(isolates, genome_id, strain_site_group)) %>%
    select(site = strain_site_group, genome_id, everything()) %>%
    filter(genome_id %in% c("g2", "g3", "g15"))



















