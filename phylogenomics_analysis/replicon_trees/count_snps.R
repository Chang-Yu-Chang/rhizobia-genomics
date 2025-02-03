#'

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv")) %>% select(contig_id, replicon_type) %>% drop_na

set_name <- "elev_med"
tt <- read_gpas(set_name)

list_sccg <- tt$gpacl %>%
    select(replicon_type, gene, genome_id) %>%
    group_by(replicon_type) %>%
    mutate(temp = T) %>%
    pivot_wider(names_from = genome_id, values_from = temp) %>%
    drop_na() %>%
    select(replicon_type, gene) %>%
    # single copy
    filter(gene %in% tt$list_sccg$gene) %>%
    ungroup()

n_genes1 <- list_sccg %>%
    group_by(replicon_type) %>%
    count(name = "n_sccg") %>%
    mutate(set_name = set_name) %>%
    relocate(set_name)


set_name <- "urbn_mel"
tt <- read_gpas(set_name)
tt$list_sccg

list_sccg <- tt$gpacl %>%
    select(replicon_type, gene, genome_id) %>%
    group_by(replicon_type) %>%
    mutate(temp = T) %>%
    pivot_wider(names_from = genome_id, values_from = temp) %>%
    drop_na() %>%
    select(replicon_type, gene) %>%
    # single copy
    filter(gene %in% tt$list_sccg$gene) %>%
    ungroup()

n_genes2 <- list_sccg %>%
    group_by(replicon_type) %>%
    count(name = "n_sccg") %>%
    mutate(set_name = set_name) %>%
    relocate(set_name)

# Number of single copy core genes
n_genes <- bind_rows(n_genes1, n_genes2)
write_csv(n_genes, paste0(folder_data, "phylogenomics_analysis/replicon_trees/n_sccg.csv"))


