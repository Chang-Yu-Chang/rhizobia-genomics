#' This script analysis the anvio generated dataset

library(tidyverse)
library(cowplot)
library(janitor)
library(FactoMineR) # for MCA
library(ggsci)
source(here::here("analysis/00-metadata.R"))

# 0. read data ----
egc <- read_delim(paste0(folder_data, "/temp/anvio/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)
#ef <- read_delim(paste0(folder_data, "/temp/anvio/enriched_functions.txt"), delim = "\t", show_col_types = F)
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

# Clean up
egc <- egc %>%
    replace_na(list(bin_name = "rest")) %>%
    filter(!is.na(bin_name)) %>%
    mutate(genome_name = factor(genome_name, c(paste0("Chang_Q5C_", 1:20), "em1021", "em1022", "wsm419"))) %>%
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
    summarize(mean(n)) # 6372


egc %>%
    #distinct(gene_cluster_id, genome_name) %>%
    filter(genome_name == "Chang_Q5C_15") %>%
    filter()

# total number of core genes
egc %>%
    distinct(gene_cluster_id, bin_name) %>%
    mutate(bin_name = factor(bin_name, c("core", "rest", "duplicate", "singleton"))) %>%
    group_by(bin_name) %>%
    count() %>%
    ungroup() %>%
    mutate(frac = n / sum(n))
# bin_name      n  frac
# 1 core       3048 0.193
# 2 rest       5808 0.367
# 3 duplicate  1896 0.120
# 4 singleton  5072 0.321




# 1. plot the Ensifer clusters using the duplicated gene pairs ----
egc %>%
    distinct(bin_name, gene_cluster_id) %>%
    filter(bin_name == "duplicate") %>%
    nrow() # 1896

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

ggsave(paste0(folder_data, "temp/44-02-gene_frequency.png"), p, width = 4, height = 4)



# 3. permutation of genome
pangenome_boots <- read_csv(paste0(folder_data, "temp/44a-pangenome_boots.csv"), show_col_types = F)
p <- pangenome_boots %>%
    pivot_longer(cols = c(-n_genomes_sampled, -bootstrap), names_to = "gene_type", values_to = "n_genes") %>%
    mutate(gene_type = factor(gene_type, c("pangene", "core", "singleton", "duplicate"))) %>%
    mutate(n_genes = n_genes / 1000) %>%
    group_by(n_genomes_sampled, gene_type) %>%
    summarize(mean_n_genes = mean(n_genes), sd_n_genes = sd(n_genes)) %>%
    ggplot() +
    geom_line(aes(x = n_genomes_sampled, y = mean_n_genes, color = gene_type), linewidth = 1) +
    geom_ribbon(aes(x = n_genomes_sampled, ymin = mean_n_genes - sd_n_genes, ymax = mean_n_genes + sd_n_genes, fill = gene_type), alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = 1:17) +
    scale_y_continuous(breaks = c(0, 5, 10, 15), minor_breaks = 0:17) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(linetype = 1),
        panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.2)
    ) +
    guides() +
    labs(x = "# of genomes included in a pangenome", y = "# of genes (k)")

ggsave(paste0(folder_data, "temp/44-03-gene_frequency_permuted.png"), plot = p, width = 6, height = 4)


# 4. permutation of genomes, scale to one
p <- pangenome_boots %>%
    mutate(across(-c(n_genomes_sampled, -bootstrap), function(x){x/pangene})) %>%  # scaled by pangene
    pivot_longer(cols = c(-n_genomes_sampled, -bootstrap), names_to = "gene_type", values_to = "n_genes") %>%
    mutate(gene_type = factor(gene_type, c("pangene", "core", "singleton", "duplicate"))) %>%
    group_by(n_genomes_sampled, gene_type) %>%
    summarize(mean_n_genes = mean(n_genes), sd_n_genes = sd(n_genes)) %>%
    ggplot() +
    geom_line(aes(x = n_genomes_sampled, y = mean_n_genes, color = gene_type), linewidth = 1) +
    geom_ribbon(aes(x = n_genomes_sampled, ymin = mean_n_genes - sd_n_genes, ymax = mean_n_genes + sd_n_genes, fill = gene_type), alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = 1:17) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0,1,0.1)) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(linetype = 1),
        panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.2)
    ) +
    guides() +
    labs(x = "# of genomes included in a pangenome", y = "# of genes (scaled)")

ggsave(paste0(folder_data, "temp/44-04-gene_frequency_permuted_scaled.png"), plot = p, width = 6, height = 4)















# 2. heatmap for genes ----
p <- gpa_long %>%
    mutate(value = factor(value)) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = gene, fill = value)) +
    scale_fill_manual(values = c(`1` = "maroon")) +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 10)
    ) +
    guides(fill = "none") +
    labs(x = "strain", y = "gene")
#ggsave(paste0(folder_data, "temp/42-02-heatmap.png"), plot = p, width = 10, height = 20)
