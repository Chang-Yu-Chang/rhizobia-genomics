#' This script analysis the anvio generated data with contig infor

library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
source(here::here("analysis/00-metadata.R"))

# 0. read data ----
# 0.1 ensifer gene presence and abense table
egc <- read_delim(paste0(folder_data, "/temp/anvio/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)

# 0.2 gene annotation
n_genomes <-  nrow(genomes)
list_annot <- rep(list(NA), n_genomes)
for (i in 1:n_genomes) {
    list_annot[[i]] <- read_delim(paste0(folder_data, "/temp/anvio/gene_annot/", genomes$genome_name[i], ".txt"), delim = "\t", show_col_types = F) %>%
        mutate(genome_id = genomes$genome_id[i], genome_name = genomes$genome_name[i])
}

gene_annot <- bind_rows(list_annot) %>%
    select(genome_name, genome_id, everything())

# 0.3 gene caller and contig
list_calls <- rep(list(NA), n_genomes)
for (i in 1:n_genomes) {
    list_calls[[i]] <- read_delim(paste0(folder_data, "/temp/anvio/gene_calls/", genomes$genome_name[i], ".txt"), delim = "\t", show_col_types = F) %>%
        mutate(genome_id = genomes$genome_id[i], genome_name = genomes$genome_name[i])
}

gene_calls <- bind_rows(list_calls) %>%
    select(genome_name, genome_id, everything())

# 0.4 create a master table with unique id for gene caller, gene cluster, genomes, and contigs
egcc <- gene_calls %>%
    left_join(egc) %>%
    mutate(genome_name = factor(genome_name, genomes$genome_name), genome_id = factor(genome_id, genomes$genome_id))

# 0.6 order the contigs by size
egcc_count <- egcc %>%
    #distinct(gene_cluster_id, gene_callers_id, .keep_all = T) %>%
    group_by(genome_name, contig) %>%
    summarize(n_genes = n()) %>%
    # Remove very small contigs
    #filter(n_genes > 10) %>%
    # Rename the contigs by size
    arrange(genome_name, desc(n_genes)) %>%
    mutate(contig_ordered = factor(paste0("c", 1:n()))) %>%
    ungroup()

# plot the
p <- egcc_count %>%
    arrange(n_genes) %>%
    ggplot() +
    geom_histogram(aes(x = n_genes), fill = "white", color = "black", binwidth = 100) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/45-00-contig_genes_n.png"), p, width = 4, height = 4)

# 0.7 label the contigs by new contig label
egcc_ordered <- egcc %>%
    left_join(select(egcc_count, genome_name, contig, contig_ordered)) %>%
    #select(genome_name, genome_id, contig, contig_ordered) %>%
    # Remove the genes on the small contigs
    filter(!is.na(contig_ordered)) %>%
    ungroup() %>%
    mutate(contig_id = paste0(genome_id, "_", contig_ordered)) %>%
    mutate(contig_id = factor(contig_id, paste0("g", rep(1:20, each = 20), "_", "c", rep(1:20, 20))))




# 1. plot the gene numbers on each contig ----
p <- egcc_count %>%
    mutate(n_genes = n_genes / 1000) %>%
    ggplot() +
    geom_col(aes(x = genome_name, y = n_genes, fill = contig_ordered), color = "black", position = position_stack(reverse = T)) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(n = 9, "Set1"), RColorBrewer::brewer.pal(n = 8, "Set2"))) +
    scale_y_continuous(breaks = seq(0, 8, 2), minor_breaks = 1:10) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(linetype = 1, linewidth = 1),
        panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.3)
    ) +
    guides(fill = guide_legend(ncol = 1, reverse = T)) +
    labs(x = "genome", y = "# of genes (k)")

ggsave(paste0(folder_data, "temp/45-01-contig_genes_n.png"), p, width = 6, height = 4)

# 2. plot the heatmap for gene presence absence
p <- egcc_ordered %>%
    filter(contig_ordered %in% paste0("c", 1:3)) %>%
    ggplot() +
    geom_rect(data = tibble(contig_ordered = paste0("c", 1:3)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = contig_ordered), alpha = 0.2) +
    geom_tile(aes(x = gene_cluster_id, y = genome_id), fill = "black") +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(n = 9, "Set1"), RColorBrewer::brewer.pal(n = 8, "Set2"))) +
    scale_y_discrete(position = "right") +
    facet_wrap(contig_ordered~., scales = "free_y", ncol = 1, strip.position = "left") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_rect(color = NA, fill = NA)
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome")


ggsave(paste0(folder_data, "temp/45-02-contig_gene.png"), p, width = 10, height = 6)

# 3. make a tree
# Compute the jaccard matrix
library(proxy)
egcc_c <- egcc_ordered %>%
    # Filter for the largest three contigs
    filter(contig_ordered %in% paste0("c", 1:3)) %>%
    select(contig_id, gene_cluster_id) %>%
    distinct(contig_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    arrange(contig_id, gene_cluster_id) %>%
    pivot_wider(id_cols = contig_id, names_from = gene_cluster_id, values_from = value, values_fill = 0)

# Compute the Jaccard distance matrix
egcc_m <- as.matrix(egcc_c[,-1])
dim(egcc_m)
rownames(egcc_m) <- as.character(egcc_c$contig_id)
jdm <- proxy::dist(egcc_m, method = "Jaccard")

# Plot the tree
library(ggtree)
te <- as.dendrogram(hclust(jdm)) # my tree!
egcc_label <- egcc_ordered %>%
    select(contig_id, genome_id, contig_ordered) %>% distinct() %>%
    filter(contig_ordered %in% paste0("c", 1:3)) %>%
    mutate(contig_id = as.character(contig_id))

p <- te %>%
    ggtree(layout = "rectangular") +
    geom_tiplab() +
    #geom_highlight(data = egcc_label, aes(node = contig_id, fill = contig_ordered)) +
    scale_x_continuous(expand = c(0.1,0)) +
    scale_y_continuous(expand = c(0.1,0)) +
    theme_tree2() +
    theme(
        plot.margin = unit(c(10, 10, 10, 10), "mm")
    ) +
    labs()


ggsave(paste0(folder_data, "temp/45-03-contig_tree.png"), p, width = 10, height = 10)


# 4. plot PCA ----
library(FactoMineR) # for MCA

# MCA because it's a categorical dataset
egcc_cm <- egcc_c %>% mutate(across(starts_with("GC"), factor))
mca <- MCA(egcc_cm[,-1], ncp = 5, graph = F)

contigs_mca <- as_tibble(mca$ind$coord) %>%
    clean_names() %>%
    bind_cols(egcc_label) %>%
    select(contig_id, everything())

mca_var <- as_tibble(mca$eig) %>%
    clean_names() %>%
    mutate(ev = eigenvalue / sum(eigenvalue))

p <- contigs_mca %>%
    ggplot() +
    geom_point(aes(x = dim_1, y = dim_2)) +
    #geom_point(aes(x = dim_1, y = dim_2, color = contig_ordered), position = position_jitter(), size = 4, stroke = 2, shape = 21, alpha = 0.8) +
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


# 5. plot the genecluster
egc %>%
    group_by(gene_cluster_id, genome_name) %>%
    count()






