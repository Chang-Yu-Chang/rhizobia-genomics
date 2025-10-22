#' This script aggregates the ani results againsts the references

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
rhizobiales <- read_delim(paste0(folder_data, "genomics/taxonomy/rhizobiales.txt"), show_col_types = F) %>%
    clean_names()

# Read ANI results
list_names <- list.files(paste0(folder_data, "genomics/taxonomy/fastani_results/")) %>%
str_subset("g\\d+_vs")

length(list_names) # should be 38 genomes * 11 batches = 418

list_files <- list()
for (i in 1:length(list_names)) {
    list_files[[i]] <- read_table(
        paste0(folder_data, "genomics/taxonomy/fastani_results/", list_names[i]),
        col_names = F, show_col_types = F
    )
    print(i)
}

ani <- bind_rows(list_files) %>%
    rename(
        query = X1,
        ref = X2,
        ani = X3,
        query_seq = X4,
        query_seq_match = X5
    ) %>%
    mutate(genome_id = str_extract(query, "g\\d+"),
           assembly_accession = str_extract(ref, "GCF_\\d+\\.\\d+")) %>%
    left_join(select(rhizobiales, organism_name, assembly_accession))


ani <- ani %>%
    select(genome_id, organism_name, ani, query_seq, query_seq_match, query, ref) %>%
    arrange(query) %>%
    group_by(query) %>%
    slice_max(ani, n = 1)

# Replace the <95% ANI with sp.
ani$organism_name[ani$ani<95] <- str_replace(ani$organism_name[ani$ani<95], " \\w+", " sp.")

write_csv(ani, paste0(folder_data, "genomics/taxonomy/ani.csv"))
