#' This script plots the pangenome

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(FactoMineR) # for MCA
library(ellipse)
library(ggsci)
source(here::here("analysis/00-metadata.R"))

#egc <- read_delim(paste0(folder_data, "genomics/pangenome/ensifer/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)
egc <- read_csv(paste0(folder_data, "temp/17-egc.csv"))
egc_wide <- read_csv(paste0(folder_data, "temp/17-egc_wide.csv"))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    filter(!genome_id %in% paste0("g", c(1,7,12,14,18)))
nrow(isolates) # This should be 32 ensifer isolates + 4 ncbi genomes

n_g <- nrow(isolates)


# 1. plot the duplicated gene pairs ----
# MCA because it's a categorical dataset
ncol(egc_wide)-1 # 14741
mca <- MCA(as.matrix(egc_wide)[,-1], ncp = 10, graph = F)

isolates_mca <- as_tibble(mca$ind$coord) %>%
    clean_names() %>%
    bind_cols(isolates) %>%
    select(genome_id, everything()) %>%
    # Manual assignment
    mutate(species = case_when(
        genome_id %in% c("usda1106", "em1021", "em1022") ~ "meliloti",
        genome_id %in% "wsm419" ~ "medicae",
        T ~ NA
    )) %>%
    mutate(rhizobia_site = ifelse(is.na(rhizobia_site), species, rhizobia_site))


mca_var <- as_tibble(mca$eig) %>%
    clean_names() %>%
    mutate(ev = eigenvalue / sum(eigenvalue))

p <- isolates_mca %>%
    ggplot() +
    geom_point(aes(x = dim_1, y = dim_2, color = rhizobia_site), size = 4, stroke = 2, shape = 21, alpha = 0.9) +
    #scale_color_d3(labels = c(H = "high elevation", L = "low elevation", ncbi = "NCBI")) +
    scale_color_manual(values = c(rhizobia_site_colors, meliloti = "grey10", medicae = "grey80"),
                       breaks = c(names(rhizobia_site_colors), "meliloti", "medicae")) +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.background = element_rect(color = NA, fill = "white"),
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides(color = guide_legend(title = "original site")) +
    labs(x = paste0("PC1 (", round(mca_var$ev[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(mca_var$ev[2]*100, 1), "%)"))

ggsave(paste0(folder_data, "temp/17a-01-mca_accessory_genes.png"), p, width = 5, height = 4)

# 2. gene frequency spectrum ----
p <- egc %>%
    group_by(num_genomes_gene_cluster_has_hits)  %>%
    distinct(gene_cluster_id, .keep_all = T) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = num_genomes_gene_cluster_has_hits, y = n), color = "black", fill = "white") +
    scale_x_continuous(breaks = seq(5, n_g, 5)) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(linetype = 2)
    ) +
    guides() +
    labs(x = "# of genomes that gene cluster has hits", y = "# of genes")

ggsave(paste0(folder_data, "temp/17a-02-gene_frequency.png"), p, width = 4, height = 4)


# 3. permutation of genome
pangenome_boots <- read_csv(paste0(folder_data, "temp/11-pangenome_boots.csv"), show_col_types = F)
p1 <- pangenome_boots %>%
    pivot_longer(cols = c(-n_genomes_sampled, -bootstrap), names_to = "bin_name", values_to = "n_genes") %>%
    mutate(bin_name = factor(bin_name, c("pangene", "core", "duplicate", "singleton", "rest"))) %>%
    mutate(n_genes = n_genes / 1000) %>%
    group_by(n_genomes_sampled, bin_name) %>%
    summarize(mean_n_genes = mean(n_genes), sd_n_genes = sd(n_genes)) %>%
    ggplot() +
    geom_line(aes(x = n_genomes_sampled, y = mean_n_genes, color = bin_name), linewidth = 1) +
    geom_ribbon(aes(x = n_genomes_sampled, ymin = mean_n_genes - sd_n_genes, ymax = mean_n_genes + sd_n_genes, fill = bin_name), alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = seq(5, n_g, 5)) +
    scale_y_continuous(breaks = seq(0,20,5), minor_breaks = 0:n_g) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(linetype = 1),
        panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.2),
        legend.position = "none"
    ) +
    guides() +
    labs(x = "# of genomes included in a pangenome", y = "# of genes (k)")

# permutation of genomes, scale to one
p2 <- pangenome_boots %>%
    mutate(across(-c(n_genomes_sampled, -bootstrap), function(x){x/pangene})) %>%  # scaled by pangene
    pivot_longer(cols = c(-n_genomes_sampled, -bootstrap), names_to = "bin_name", values_to = "n_genes") %>%
    mutate(bin_name = factor(bin_name, c("pangene", "core", "duplicate", "singleton", "rest"))) %>%
    group_by(n_genomes_sampled, bin_name) %>%
    summarize(mean_n_genes = mean(n_genes), sd_n_genes = sd(n_genes)) %>%
    ggplot() +
    geom_line(aes(x = n_genomes_sampled, y = mean_n_genes, color = bin_name), linewidth = 1) +
    geom_ribbon(aes(x = n_genomes_sampled, ymin = mean_n_genes - sd_n_genes, ymax = mean_n_genes + sd_n_genes, fill = bin_name), alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = seq(5, n_g, 5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0,1,0.1)) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(linetype = 1),
        panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.2)
    ) +
    guides() +
    labs(x = "# of genomes included in a pangenome", y = "# of genes (scaled)")

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 1.2), labels = c("A", "B"), scale = 0.95) + paint_white_background()

ggsave(paste0(folder_data, "temp/17a-03-gene_frequency_permuted.png"), plot = p, width = 12, height = 4)

#
egc_count <- egc %>%
    distinct(bin_name, gene_cluster_id) %>%
    mutate(bin_name = factor(bin_name, c("pangene", "core", "duplicate", "singleton", "rest"))) %>%
    group_by(bin_name) %>%
    count()

egc_count %>%
    ungroup() %>%
    mutate(frac = n / sum(n))

sum(egc_count$n) # 17729





















