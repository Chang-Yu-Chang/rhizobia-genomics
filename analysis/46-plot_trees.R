#' This script plot the trees for allthe 17 genomes based on many measures

library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
# library(FactoMineR) # for MCA
library(proxy) # for computing jaccard matrix
library(ape)
library(ggtree)
source(here::here("analysis/00-metadata.R"))

# 1. tree based on fastani ----
tb_ani <- read_csv(paste0(folder_data, "temp/41-tb_ani.csv"), show_col_types = F)
tb_ani <- tb_ani %>%
    left_join(rename(genomes, g_A = genome_name, id_A = genome_id)) %>%
    left_join(rename(genomes, g_B = genome_name, id_B = genome_id)) %>%
    select(id_A, id_B, ani)

tb_ani_m <- tb_ani %>%
    pivot_wider(id_cols = id_B, names_from = id_A, values_from = ani)
tb_ani_m <- as.dist(tb_ani_m[,-1])
te <- as.dendrogram(hclust(tb_ani_m))

p <- te %>%
    ggtree(layout = "rectangular") +
    geom_tiplab() +
    scale_x_continuous(expand = c(0.1,0)) +
    scale_y_continuous(expand = c(0.1,0)) +
    theme_tree2() +
    theme(
        plot.margin = unit(c(10, 10, 10, 10), "mm")
    ) +
    labs()
ggsave(paste0(folder_data, "temp/46-01-genome_ani.png"), p, width = 10, height = 10)

# 2. tree based on whole genome gene presence an absence ----
egcct <- read_csv(paste0(folder_data, "temp/45-egcct.csv"), show_col_types = F)
egcct_c <- egcct %>%
    select(genome_id, gene_cluster_id) %>%
    distinct(genome_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    arrange(genome_id, gene_cluster_id) %>%
    pivot_wider(id_cols = genome_id, names_from = gene_cluster_id, values_from = value, values_fill = 0)

egcc_m <- as.matrix(egcct_c[,-1])
dim(egcc_m)
rownames(egcc_m) <- as.character(egcct_c$genome_id)
jdm <- proxy::dist(egcc_m, method = "Jaccard")

# Plot the tree
te <- as.dendrogram(hclust(jdm))
p <- te %>%
    ggtree(layout = "rectangular") +
    geom_tiplab() +
    scale_x_continuous(expand = c(0.1,0)) +
    scale_y_continuous(expand = c(0.1,0)) +
    theme_tree2() +
    theme(
        plot.margin = unit(c(10, 10, 10, 10), "mm")
    ) +
    labs()

ggsave(paste0(folder_data, "temp/46-02-genome_pa.png"), p, width = 10, height = 10)

# 3. tree based on chromosomes ---
egcct <- read_csv(paste0(folder_data, "temp/45-egcct.csv"), show_col_types = F)
compute_jaccard <- function (egcct, contig_t) {
    egcct_c <- egcct %>%
        filter(contig_type == contig_t) %>%
        select(contig_id, gene_cluster_id) %>%
        distinct(contig_id, gene_cluster_id) %>%
        mutate(value = 1) %>%
        arrange(contig_id, gene_cluster_id) %>%
        pivot_wider(id_cols = contig_id, names_from = gene_cluster_id, values_from = value, values_fill = 0)

    egcc_m <- as.matrix(egcct_c[,-1])
    dim(egcc_m)
    rownames(egcc_m) <- as.character(egcct_c$contig_id)
    jdm <- proxy::dist(egcc_m, method = "Jaccard")
    return(jdm)
}
make_tree <- function(jdm) as.phylo(hclust(jdm))
plot_tree <- function(tree, contig_t) {
    tree %>%
        ggtree(layout = "rectangular") +
        geom_tiplab() +
        scale_x_continuous(expand = c(0.1,0)) +
        scale_y_continuous(expand = c(0.1,0)) +
        theme_tree2() +
        theme(
            plot.margin = unit(c(10, 10, 10, 10), "mm")
        ) +
        labs(title = contig_t)
}

tree_chr <- compute_jaccard(egcct, "chromosome") %>% make_tree()
p <- plot_tree(tree_chr, "chromosome")
ggsave(paste0(folder_data, "temp/46-03-chromosome_pa.png"), p, width = 10, height = 10)

# 4. tree based on psymA ----
tree_a <- compute_jaccard(egcct, "psyma") %>% make_tree()
p <- plot_tree(tree_a, "pSymA")
ggsave(paste0(folder_data, "temp/46-04-psyma_pa.png"), p, width = 10, height = 10)

# 5. tree based on psymB ----
tree_b <- compute_jaccard(egcct, "psymb") %>% make_tree()
p <- plot_tree(tree_b, "pSymB")
ggsave(paste0(folder_data, "temp/46-05-psymb_pa.png"), p, width = 10, height = 10)




# Compute the Colless index
library(treebalance) # for calculating Colless index
collessI(tree_chr)
collessI(tree_a)
collessI(tree_b)






























