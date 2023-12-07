#' This script analyzes the pangenome data

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    source(here::here("analysis/00-metadata.R"))
})

egc <- read_delim(paste0(folder_data, "genomics/pangenome/ensifer/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    filter(!genome_id %in% paste0("g", c(1,7,12,14,18)))
nrow(isolates) # This should be 32 ensifer isolates + 4 ncbi genomes


# 0. Clean up ----
egc <- clean_names(egc)

# Count the gene cluster by the number of genomes they appear
n_total_genomes <- length(unique(egc$genome_name))
nrow(n_total_genomes) # 36
egcc <- egc %>%
    distinct(gene_cluster_id, genome_name) %>%
    group_by(gene_cluster_id) %>%
    count(name = "n_genomes") %>%
    ungroup %>%
    mutate(bin_name = case_when(
        n_genomes == n_total_genomes ~ "core",
        n_genomes == 1 ~ "singleton",
        n_genomes == 2 ~ "duplicate",
        T ~ "rest"
    ))

egc <- egc %>%
    select(-bin_name) %>%
    left_join(egcc) %>%
    mutate(genome_name = factor(genome_name, isolates$genome_name)) %>%
    arrange(genome_name) %>%
    mutate(bin_name = factor(bin_name, c("core", "rest", "duplicate", "singleton"))) %>%
    left_join(isolates) %>%
    select(unique_id, gene_cluster_id, bin_name, genome_name, genome_id, everything())


write_csv(egc, paste0(folder_data, "temp/17-egc.csv"))

# 1. numbers ----
# Total number of gene clusters
egc %>%
    distinct(gene_cluster_id) %>%
    nrow() # 17729
# mean of total number of gene cluster in a genome
egc %>%
    distinct(gene_cluster_id, genome_name) %>%
    group_by(genome_name) %>%
    count() %>%
    ungroup() %>%
    summarize(mean(n)) # 6308


# total number of core genes
egc %>%
    group_by(bin_name) %>%
    count() %>%
    ungroup() %>%
    mutate(frac = n / sum(n))
# bin_name       n   frac
# <fct>      <int>  <dbl>
# 1 core      121929 0.465
# 2 rest      131324 0.501
# 3 duplicate   3806 0.0145
# 4 singleton   5286 0.0201


# 1. clustering using accessory genes ----
egc %>%
    distinct(bin_name, gene_cluster_id) %>%
    filter(bin_name != "core") %>%
    nrow() # 14741

# Reshape
egc_wide <- egc %>%
    # Remove the core genes because they are not helpful
    filter(bin_name != "core") %>%
    distinct(gene_cluster_id, genome_name, .keep_all = T) %>%
    select(gene_cluster_id, genome_name) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = gene_cluster_id, values_from = value, values_fill = 0) %>%
    mutate(across(starts_with("GC"), factor))

write_csv(egc_wide, paste0(folder_data, "temp/17-egc_wide.csv"))

# 2. permutation ----
# Extract only the unique gene id
egc_g <- egc %>%
    select(unique_id, gene_cluster_id, genome_name) %>%
    distinct(gene_cluster_id, genome_name, .keep_all = T) %>%
    left_join(isolates) %>%
    select(-genome_name)

n_b = 100 # number of permutation
n_g = 36 # total number of genomes in the pangenome

list_pangenome <- rep(list(NA), n_g-1)

for (j in 1:(n_g-1)) { # number of genomes sampled
    list_boot <- rep(list(NA), n_b)
    for (i in 1:n_b) {
        set.seed(i)
        g_sampled <- sample(isolates$genome_id, size = j, replace = F)

        gfs <- egc_g %>% # gene frequency spectrum
            filter(genome_id %in% g_sampled) %>%
            group_by(gene_cluster_id) %>%
            summarize(n_genomes = n()) %>%
            group_by(n_genomes) %>%
            summarize(n_genes = n())

        g1 <- gfs$n_genes[gfs$n_genomes == 1] # singleton
        g2 <- gfs$n_genes[gfs$n_genomes == 2] # duplicate
        gc <- gfs$n_genes[gfs$n_genomes == j] # core gene clusters
        gp <- sum(gfs$n_genes) # total number of gene clusters in the pangenome

        list_boot[[i]] = tibble(bootstrap = i, singleton = g1, duplicate = g2, core = gc, pangene = gp)
    }

    list_pangenome[[j]]  <- bind_rows(list_boot) %>%
        mutate(n_genomes_sampled = j) %>%
        select(n_genomes_sampled, everything())

    print(j)
}

pangenome_boots <- bind_rows(list_pangenome)

write_csv(pangenome_boots, file = paste0(folder_data, "temp/17-pangenome_boots.csv"))













