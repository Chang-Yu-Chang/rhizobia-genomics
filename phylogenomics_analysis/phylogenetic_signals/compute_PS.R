#' This scripts implements the phylogenetic comparative analysis

renv::load()
library(tidyverse)
library(phytools)
source(here::here("metadata.R"))

# Traits
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_traits <- read_csv(paste0(folder_data, "phenotypes_analysis/isolates_traits.csv"))
isolates <- isolates %>%
    left_join(isolates_contigs) %>%
    left_join(isolates_traits) %>%
    filter(!genome_id %in% c("g28"))

list_traits <- c(paste0(rep(c("r", "lag", "maxOD"), each = 4), "_", rep(c(25, 30, 35, 40), 3), "c"),
                 "root_biomass_mg", "shoot_biomass_mg", "nodule_count")

# Trees
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

# Compute phylogenetic signals
list_traits <- str_subset(names(isolates), "0c$|5c$|nodule|biomass|contig_length")
compute_ps <- function(tree, genome_id, trait_value) {
    set.seed(1)
    temp <- setNames(trait_value, genome_id)
    ps1 <- phylosig(tree, temp, test = T, nsim = 1000, method = "K")
    ps2 <- phylosig(tree, temp, test = T, nsim = 1000, method = "lambda")
    return(tibble(k = ps1$K, p_k = ps1$P, lambda = ps2$lambda, p_lambda = ps2$P))
}

tb1 <- tibble(trait = list_traits) %>%
    rowwise() %>%
    mutate(result = list(compute_ps(tr, isolates$genome_id, isolates[[trait]]))) %>%
    unnest(result)

tb2 <- tibble(trait = list_traits) %>%
    rowwise() %>%
    mutate(result = list(compute_ps(tr_gpa, isolates$genome_id, isolates[[trait]]))) %>%
    unnest(result)

traits_ps <- bind_rows(mutate(tb1, tree = "core"), mutate(tb2, tree = "gpa"))
write_csv(traits_ps, paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))
