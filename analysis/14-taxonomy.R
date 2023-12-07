#' This script assigns the taxonomy to isolates

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F)

# 1. aggregate mash results ----
# 1.1 mash on whole genome ----
list_screen <- rep(list(NA), nrow(isolates))
for (i in 1:nrow(isolates)) {
    mash_dist <- read_tsv(paste0(folder_genomes, isolates$genome_name[i], "/04-taxonomy/mash/screen.tab"), show_col_types = F, col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment"))
    list_screen[[i]] <- mash_dist %>%
        arrange(p_value) %>%
        mutate(genome_id = isolates$genome_id[i]) %>%
        select(genome_id, everything())
}

mash_g <- list_screen[-1] %>% bind_rows()

# Clean up
mash_g <- mash_g %>%
    left_join(isolates) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))

write_csv(mash_g, paste0(folder_data, "temp/14-mash_g.csv"))

# 1.2 mash on individual contigs
list_screen <- rep(list(NA), nrow(isolates))

for (i in 1:nrow(isolates)) {
    list_contig_names <- list.files(paste0(folder_genomes, isolates$genome_name[i], "/04-taxonomy/mash/"), pattern = "c_") %>%
        str_subset(".tab") %>% str_remove(".tab") %>% unique()
    list_contig_screen <- rep(list(NA), length(list_contig_names))
    for (j in 1:length(list_contig_names)) {
        mash_dist <- read_tsv(paste0(folder_genomes, isolates$genome_name[i], "/04-taxonomy/mash/", list_contig_names[j], ".tab"), show_col_types = F, col_names = c("identity", "shared_hashes", "median_multiplicity", "p_value", "query_ID", "query_comment"))
        if (nrow(mash_dist)==0) next
        list_contig_screen[[j]] <- mash_dist %>%
            arrange(p_value) %>%
            mutate(contig_id = list_contig_names[j]) %>%
            select(contig_id, everything())
    }
    list_contig_screen <- list_contig_screen[!is.na(list_contig_screen)]
    list_screen[[i]] <- bind_rows(list_contig_screen) %>%
        mutate(genome_id = isolates$genome_id[i]) %>%
        select(genome_id, contig_id, everything())
}
mash_c <- list_screen %>% bind_rows()
write_csv(mash_c, paste0(folder_data, "temp/14-mash_c.csv"))

# 2. filter the mash top hits ----
# 2.1 genome top hits ----
mash_g_top <- mash_g %>%
    group_by(genome_id) %>%
    slice(1:3)
write_csv(mash_g_top, paste0(folder_data, "temp/14-mash_g_top.csv"))


# 2.2 contig top hits ----
mash_c_top <- mash_c %>%
    group_by(genome_id, contig_id) %>%
    slice(1:3)
write_csv(mash_c_top, paste0(folder_data, "temp/14-mash_c_top.csv"))


# 3. clean the taxonomy name ----
# Genomes
isolates_mash <- mash_g_top %>%
    slice(1) %>%
    select(genome_id, genome_name, exp_id, rhizobia_site, query_comment) %>%
    # Extract the species name
    mutate(species_name = str_extract(query_comment, "Sinorhizobium.*|Ensifer.*|Rhizobium.*")) %>%
    mutate(species_name = str_extract(species_name, "\\S+\\s+\\S+")) %>%
    mutate(species_name = str_replace(species_name, "Sinorhizobium", "Ensifer")) %>%
    # Genus0
    mutate(genus = str_extract(species_name, "\\S+"))

write_csv(isolates_mash, paste0(folder_data, "temp/14-isolates_mash.csv"))

# Contigs
contigs_large <- read_csv(paste0(folder_data, "temp/12-contigs_large.csv"), show_col_types = F)
contigs_mash <- mash_c_top %>%
    slice(1) %>%
    right_join(contigs_large) %>%
    left_join(isolates) %>%
    select(genome_id, contig_length, query_comment) %>%
    # Assign genetic element by size
    group_by(genome_id) %>%
    mutate(ge_type = case_when(
        contig_length == max(contig_length) ~ "chromosome",
        T ~ "plasmid"
    )) %>%
    # Genetic element name
    mutate(ge_name = ifelse(ge_type == "plasmid", str_extract(query_comment, ".*plasmid\\s\\S+") %>% str_remove(".*plasmid\\s") %>% str_remove(","), NA)) %>%
    mutate(ge_name = ifelse(ge_type == "chromosome", "chromosome", ge_name)) %>%
    # # Extract the contig name
    # mutate(contig_name = str_extract(query_comment, "Ensifer\\s\\S+|Sinorhizobium.*plasmid\\s\\S+|Sinorhizobium.*chromosome|Sinorhizobium\\s\\S+\\s\\S+|Rhizobium.*plasmid|Rhizobium\\s\\S+\\s\\S+")) %>%
    # # Extract genus
    # mutate(genus = str_extract(contig_name, "\\S+") %>% str_replace("Sinorhizobium", "Ensifer")) %>%
    # # Extract species name
    # mutate(species = str_extract(contig_name, "\\S+\\s\\S+") %>% str_replace("Sinorhizobium", "Ensifer")) %>%
    # Genetic element
    select(genome_id, contig_id, contig_length,  ge_type, ge_name) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id, desc(contig_length))


write_csv(contigs_mash, paste0(folder_data, "temp/14-contigs_mash.csv"))

# # 3. find the contigs with top length
# mash_c_top <- read_csv(paste0(folder_data, "temp/38-mash_c_top.csv"), show_col_types = F)
#
# egcct <- read_csv(paste0(folder_data, "temp/45-egcct.csv"), show_col_types = F)
# egcct_c <- egcct %>%
#     distinct(genome_id, contig, contig_ordered, contig_type)
#
# mash_c_top %>%
#     select(genome_id, contig = contig_id, identity, query_comment) %>%
#     left_join(egcct_c) %>%
#     mutate(genome_id = factor(genome_id, paste0("g", 1:20))) %>%
#     arrange(genome_id, contig_type) %>%
#     drop_na() %>%
#     #filter(str_detect(query_comment, "Ensifer|Sinorhizobium")) %>%
#     select(genome_id, contig_type, query_comment) %>%
#     group_by(genome_id, contig_type) %>%
#     mutate(mash_hits = paste0("mash", 1:n())) %>%
#     ungroup() %>%
#     pivot_wider(id_cols = c(genome_id, contig_type), names_from = mash_hits, values_from = query_comment) %>%
#     left_join(select(isolates, genome_id, strain_site_group)) %>%
#     select(site = strain_site_group, genome_id, everything()) %>%
#     filter(genome_id %in% c("g2", "g3", "g15"))
# #
#
# # 3. aggregate sourmash results ----
# # 3.1 read sourmahs results
# list_sm <- rep(list(NA), nrow(isolates))
#
# for (i in 1:length(list_sm)) {
#     list_sm[[i]] <- read_csv(paste0(folder_genomes, isolates$genome_name, "/04-taxonomy/sourmash/gathered.csv"), show_col_types = F) %>%
#         mutate(query_md5 = as.character(query_md5)) %>%
#         mutate(genome_id = i)
# }
#
# tb_sm <- bind_rows(list_sm[-1]) %>%  select(genome_id, everything())
# write_csv(tb_sm, paste0(folder_data, "temp/14-tb_sm.csv"))
#
# tb_sm %>%
#     group_by(genome_id) %>%
#     summarize(n_matches = n())
# # 3.2
# tb_sm_tax <- tb_sm %>%
#     group_by(genome_id) %>%
#     slice(1) %>%
#     #select(genome_id, intersect_bp, name) %>%
#     mutate(`intersect_bp (Mbp)` = intersect_bp/1000000) %>%
#     mutate(name = str_replace(name, "GC[A|F]_\\d+\\.1 ", "")) %>%
#     #mutate(contig_id = 1:n()) %>%
#     select(genome_id, `intersect_bp (Mbp)`, f_orig_query, f_match, name)
#
# write_csv(tb_sm_tax, paste0(folder_data, "temp/14-tb_sm_tax.csv"))
#
#
# # 4. extract the genus and species level information from mash ----



