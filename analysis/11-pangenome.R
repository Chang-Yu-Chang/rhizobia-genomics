#' This script analyzes the pangenome data

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    library(FactoMineR) # for MCA
    library(ggsci)
    source(here::here("analysis/00-metadata.R"))
})

egc <- read_delim(paste0(folder_data, "genomics/pangenome/ensifer/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)
genomes_mapping <- read_csv(paste0(folder_data, "temp/00-genomes_mapping.csv"), show_col_types = F)
isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"), show_col_types = F)
isolates <- full_join(genomes_mapping, isolates_mapping) %>%
    select(genome_name, genome_id, exp_id, rhizobia_site) %>%
    filter(!genome_id %in% paste0("g", c(1,7,12,14,18)))
nrow(isolates) # This should be 32 + 4


# Clean up
egc <- clean_names(egc)

## Count the number of genomes that each gene cluster appears
n_total_genomes <- length(unique(egc$genome_name))
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
    arrange(genome_name)


# Numbers
# Total number of gene clusters
egc %>%
    distinct(gene_cluster_id) %>%
    nrow() # 15824
# mean of total number of gene cluster in a genome
egc %>%
    distinct(gene_cluster_id, genome_name) %>%
    group_by(genome_name) %>%
    count() %>%
    ungroup() %>%
    summarize(mean(n)) # 6337


# total number of core genes
egc %>%
    distinct(gene_cluster_id, bin_name) %>%
    mutate(bin_name = factor(bin_name, c("core", "rest", "duplicate", "singleton"))) %>%
    group_by(bin_name) %>%
    count() %>%
    ungroup() %>%
    mutate(frac = n / sum(n))
# bin_name      n  frac
# <fct>     <int> <dbl>
# 1 core       3048 0.193
# 2 rest       6054 0.383
# 3 duplicate  1794 0.113
# 4 singleton  4928 0.311




# 1. plot the Ensifer clusters using the duplicated gene pairs ----
egc %>%
    distinct(bin_name, gene_cluster_id) %>%
    filter(bin_name == "duplicate") %>%
    nrow() # 1794

# Reshape
egc_wide <- egc %>%
    # Remove the core genes because they are not helpful
    filter(bin_name == "duplicate") %>%
    distinct(gene_cluster_id, genome_name, .keep_all = T) %>%
    select(gene_cluster_id, genome_name) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = gene_cluster_id, values_from = value, values_fill = 0) %>%
    mutate(across(starts_with("GC"), factor))


# MCA because it's a categorical dataset
mca <- MCA(egc_wide[,-1], ncp = 5, graph = F)

isolates_mca <- as_tibble(mca$ind$coord) %>%
    clean_names() %>%
    bind_cols(isolates[c(1:18),]) %>%
    select(genome_id, everything())

mca_var <- as_tibble(mca$eig) %>%
    clean_names() %>%
    mutate(ev = eigenvalue / sum(eigenvalue))

p <- isolates_mca %>%
    ggplot() +
    geom_point(aes(x = dim_1, y = dim_2, color = rhizobia_site), size = 4, stroke = 2, shape = 21, alpha = 0.8) +
    scale_color_d3(labels = c(H = "high elevation", L = "low elevation", ncbi = "NCBI")) +
    theme_classic() +
    theme(
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(color = "black", fill = "white")
    ) +
    guides(color = guide_legend(title = "original site")) +
    labs(x = paste0("PC1 (", round(mca_var$ev[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(mca_var$ev[2]*100, 1), "%)"))

ggsave(paste0(folder_data, "temp/11-01-duplicated_gene_pairs.png"), p, width = 4, height = 4)


# 2. gene frequency spectrum ----
p <- egc %>%
    group_by(num_genomes_gene_cluster_has_hits)  %>%
    distinct(gene_cluster_id, .keep_all = T) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = num_genomes_gene_cluster_has_hits, y = n), color = "black", fill = "white") +
    scale_x_continuous(breaks = 1:17) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(linetype = 2)
    ) +
    guides() +
    labs(x = "# of genomes that gene cluster has hits", y = "# of genes")

ggsave(paste0(folder_data, "temp/11-02-gene_frequency.png"), p, width = 4, height = 4)

#
#
# # 3. permutation of genome
# pangenome_boots <- read_csv(paste0(folder_data, "temp/11a-pangenome_boots.csv"), show_col_types = F)
# p <- pangenome_boots %>%
#     pivot_longer(cols = c(-n_genomes_sampled, -bootstrap), names_to = "gene_type", values_to = "n_genes") %>%
#     mutate(gene_type = factor(gene_type, c("pangene", "core", "singleton", "duplicate"))) %>%
#     mutate(n_genes = n_genes / 1000) %>%
#     group_by(n_genomes_sampled, gene_type) %>%
#     summarize(mean_n_genes = mean(n_genes), sd_n_genes = sd(n_genes)) %>%
#     ggplot() +
#     geom_line(aes(x = n_genomes_sampled, y = mean_n_genes, color = gene_type), linewidth = 1) +
#     geom_ribbon(aes(x = n_genomes_sampled, ymin = mean_n_genes - sd_n_genes, ymax = mean_n_genes + sd_n_genes, fill = gene_type), alpha = 0.2, color = NA) +
#     scale_x_continuous(breaks = 1:17) +
#     scale_y_continuous(breaks = c(0, 5, 10, 15), minor_breaks = 0:17) +
#     theme_classic() +
#     theme(
#         panel.grid.major.y = element_line(linetype = 1),
#         panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.2)
#     ) +
#     guides() +
#     labs(x = "# of genomes included in a pangenome", y = "# of genes (k)")
#
# ggsave(paste0(folder_data, "temp/11-03-gene_frequency_permuted.png"), plot = p, width = 6, height = 4)
#
#
# # 4. permutation of genomes, scale to one
# p <- pangenome_boots %>%
#     mutate(across(-c(n_genomes_sampled, -bootstrap), function(x){x/pangene})) %>%  # scaled by pangene
#     pivot_longer(cols = c(-n_genomes_sampled, -bootstrap), names_to = "gene_type", values_to = "n_genes") %>%
#     mutate(gene_type = factor(gene_type, c("pangene", "core", "singleton", "duplicate"))) %>%
#     group_by(n_genomes_sampled, gene_type) %>%
#     summarize(mean_n_genes = mean(n_genes), sd_n_genes = sd(n_genes)) %>%
#     ggplot() +
#     geom_line(aes(x = n_genomes_sampled, y = mean_n_genes, color = gene_type), linewidth = 1) +
#     geom_ribbon(aes(x = n_genomes_sampled, ymin = mean_n_genes - sd_n_genes, ymax = mean_n_genes + sd_n_genes, fill = gene_type), alpha = 0.2, color = NA) +
#     scale_x_continuous(breaks = 1:17) +
#     scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0,1,0.1)) +
#     theme_classic() +
#     theme(
#         panel.grid.major.y = element_line(linetype = 1),
#         panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.2)
#     ) +
#     guides() +
#     labs(x = "# of genomes included in a pangenome", y = "# of genes (scaled)")
#
# ggsave(paste0(folder_data, "temp/11-04-gene_frequency_permuted_scaled.png"), plot = p, width = 6, height = 4)













