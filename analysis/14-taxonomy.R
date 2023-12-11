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

# 4. aggregate 16S blast results ----
# 4.1 blast on all 16S ----
list_blasts <- rep(list(NA), nrow(isolates))
for (i in 1:nrow(isolates)) {
    blast <- read_delim(paste0(folder_genomes, isolates$genome_name[i], "/04-taxonomy/16s/blast.txt"), show_col_types = F, delim = "\t", col_names = c("qseqid", "sseqid", "stitle", "bitscore", "evalue", "length", "pident"))
    list_blasts[[i]] <- blast %>%
        arrange(qseqid, desc(bitscore)) %>%
        mutate(genome_id = isolates$genome_id[i]) %>%
        select(genome_id, everything())
}

blast_16s <- bind_rows(list_blasts)

# Clean up
blast_16s <- blast_16s %>%
    mutate(species_name = str_extract(stitle, "\\S+\\s\\S+\\s\\S+") %>% str_remove("NR_\\S+\\s")) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    group_by(genome_id, qseqid) %>%
    # πercentage of identical matches >99%
    filter(pident > 99) %>%
    distinct(species_name, .keep_all = T) %>%
    left_join(isolates)

write_csv(blast_16s, paste0(folder_data, "temp/14-blast_16s.csv"))
#
# blast_16s %>%
#     filter(genome_id == "g8") %>%
#     select(genome_id, qseqid, stitle) %>%
#     view
