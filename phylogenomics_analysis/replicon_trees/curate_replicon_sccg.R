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

list_sccg_chromosome <- list_sccg %>% filter(replicon_type == "chromosome") %>% select(gene)
list_sccg_pSymA <- list_sccg %>% filter(replicon_type == "pSymA") %>% select(gene)
list_sccg_pSymB <- list_sccg %>% filter(replicon_type == "pSymB") %>% select(gene)
list_sccg_pAcce <- list_sccg %>% filter(replicon_type == "pAcce") %>% select(gene) # 8 core on pAcce

write_csv(list_sccg, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg.csv"))
write_csv(list_sccg_chromosome, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg_chromosome.csv"), col_names = F)
write_csv(list_sccg_pSymA, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg_pSymA.csv"), col_names = F)
write_csv(list_sccg_pSymB, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg_pSymB.csv"), col_names = F)
write_csv(list_sccg_pAcce, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg_pAcce.csv"), col_names = F)

n_genes1 <- list_sccg %>%
    group_by(replicon_type) %>%
    count(name = "n_sccg") %>%
    mutate(set_name = set_name) %>%
    relocate(set_name)

#
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

list_sccg_chromosome <- list_sccg %>% filter(replicon_type == "chromosome") %>% select(gene)
list_sccg_pSymA <- list_sccg %>% filter(replicon_type == "pSymA") %>% select(gene)
list_sccg_pSymB <- list_sccg %>% filter(replicon_type == "pSymB") %>% select(gene)

write_csv(list_sccg, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg.csv"))
write_csv(list_sccg_chromosome, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg_chromosome.csv"))
write_csv(list_sccg_pSymA, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg_pSymA.csv"), col_names = F)
write_csv(list_sccg_pSymB, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/list_sccg_pSymB.csv"), col_names = F)

n_genes2 <- list_sccg %>%
    group_by(replicon_type) %>%
    count(name = "n_sccg") %>%
    mutate(set_name = set_name) %>%
    relocate(set_name)

# Num ber of single copy core genes
n_genes <- bind_rows(n_genes1, n_genes2)
write_csv(n_genes, paste0(folder_data, "phylogenomics_analysis/replicon_trees/n_sccg.csv"))


