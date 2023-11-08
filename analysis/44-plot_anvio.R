#' This script analysis the anvio generated dataset

library(tidyverse)
library(cowplot)
library(janitor)
library(FactoMineR) # for MCA
library(ggsci)
source(here::here("analysis/00-metadata.R"))

# 0. read data
egc <- read_delim(paste0(folder_data, "/temp/anvio/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)
ef <- read_delim(paste0(folder_data, "/temp/anvio/enriched_functions.txt"), delim = "\t", show_col_types = F)
isolates_anvio <- read_delim(paste0(folder_data, "temp/42-isolates_anvio.txt"), show_col_types = F) %>%
    select(sample, r_category)
isolates <- read_csv(paste0(folder_data, "temp/42-isolates.csv"), show_col_types = F) %>%
    mutate(sample = str_replace(genome_id, "g", "Chang_Q5C_")) %>%
    select(sample, everything(), -r, -lag, -maxOD) %>%
    bind_rows(tibble(sample = c("em1021", "em1022", "wsm419"))) %>%
    left_join(isolates_anvio) %>%
    mutate(strain = ifelse(is.na(strain), "ncbi", strain),
           strain_site = ifelse(is.na(strain_site), "ncbi", strain_site),
           strain_site_group = ifelse(is.na(strain_site_group), "ncbi", strain_site_group))

# Numbers
# Total number of gene clusters
egc %>%
    distinct(gene_cluster_id) %>%
    nrow() # 10752
# mean of total number of gene cluster in a genome
egc %>%
    distinct(gene_cluster_id, genome_name) %>%
    group_by(genome_name) %>%
    count() %>%
    ungroup() %>%
    summarize(mean(n)) # 6074
# total number of core genes
egc %>%
    distinct(gene_cluster_id, bin_name) %>%
    replace_na(list(bin_name = "rest")) %>%
    group_by(bin_name) %>%
    count()
# bin_name                 n
# <chr>                <int>
# 1 better_core            654
# 2 core                  1658
# 3 duplicated_gene_pair  1721
# 4 rest                  6719



# Clean up
egc <- egc %>%
    filter(!is.na(bin_name)) %>%
    mutate(genome_name = factor(genome_name, c(paste0("Chang_Q5C_", 1:20), "em1021", "em1022", "wsm419"))) %>%
    arrange(genome_name)

egc_wide <- egc %>%
    # Remove the core genes because they are not helpful
    filter(bin_name == "duplicated_gene_pair") %>%
    select(gene_cluster_id, genome_name) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = gene_cluster_id, values_from = value, values_fill = 0, ) %>%
    mutate(across(starts_with("GC"), factor))

# Numbers
egc %>%
    distinct(gene_cluster_id, genome_name) %>%
    filter()


egc %>%
    mutate(bin_name = factor(bin_name, c("core", "better_core", "duplicated_gene_pair"))) %>%
    group_by(bin_name) %>%
    summarize(count = n())

# 1. plot the Ensifer clusters using the duplicated gene pairs
egc %>%
    filter(bin_name == "duplicated_gene_pair") %>%
    nrow()
# 3442

# MCA because it's a categorical dataset
mca <- MCA(egc_wide[,-1], ncp = 5, graph = F)

isolates_mca <- as_tibble(mca$ind$coord) %>%
    clean_names() %>%
    bind_cols(isolates) %>%
    select(sample, everything())

mca_var <- as_tibble(mca$eig) %>%
    clean_names() %>%
    mutate(ev = eigenvalue / sum(eigenvalue))

p <- isolates_mca %>%
    ggplot() +
    geom_point(aes(x = dim_1, y = dim_2, color = strain_site_group), size = 4, stroke = 2, shape = 21, alpha = 0.8) +
    scale_color_d3(labels = c(H = "high elevation", L = "low elevation", ncbi = "NCBI")) +
    theme_classic() +
    theme(
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(color = "black", fill = "white")
    ) +
    guides(color = guide_legend(title = "original site")) +
    labs(x = paste0("PC1 (", round(mca_var$ev[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(mca_var$ev[2]*100, 1), "%)"))

ggsave(paste0(folder_data, "temp/44-01-duplicated_gene_pairs.png"), p, width = 4, height = 4)

# 2. Plot the enriched functions
ef %>%
    filter(unadjusted_p_value < 0.05) %>%
    ggplot() +
    geom_histogram(aes(x = enrichment_score)) +
    theme_classic() +
    theme() +
    guides() +
    labs()
