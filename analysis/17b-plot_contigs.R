#' This script plots the gene content by the contigs

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    library(ggsci)
    library(proxy) # for computing jaccard matrix
    library(ape)
    library(ggtree)
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    filter(!genome_id %in% paste0("g", c(1,7,12,14,18)))
egc <- read_csv(paste0(folder_data, "temp/17-egc.csv"), show_col_types = F)
gcalls <- read_csv(paste0(folder_data, "temp/17-gcalls.csv"), show_col_types = F)
contigs_mash <- read_csv(paste0(folder_data, "temp/14-contigs_mash.csv"), show_col_types = F)
isolates_mash <- read_csv(paste0(folder_data, "temp/14-isolates_mash.csv"), show_col_types = F)
nrow(isolates) # This should be 32 ensifer isolates + 4 ncbi genomes

#n_g <- nrow(isolates)

# 0. Data wrangling ----
# 0.1.
contigs_mash <- contigs_mash %>%
    drop_na() %>%
    unite(col = "contig_unique_id", genome_id, contig_id, remove = F)


# 0.2. join data ----
egcalls <- egc %>%
    left_join(gcalls) %>%
    select(unique_id, gene_cluster_id, genome_id, contig_id, gene_callers_id, max_num_paralogs, prokka_prodigal_acc) %>%
    left_join(contigs_mash) %>%
    right_join(isolates) %>%
    filter(!is.na(contig_length))


# 0.3 gene numbers on each ----
# egcalls %>%
#     unite(col = "contig_unique_id", genome_id, contig_id, remove = F) %>%
#     select(gene_cluster_id, contig_unique_id) %>%
#     distinct() %>%
#     mutate(value = 1) %>%
#     pivot_wider(names_from = contig_unique_id, values_fill = 0) %>%
#     filter()


# 1. plot the number of genes on each contig ----
p <- egcalls %>%
    distinct(genome_id, contig_id, contig_length, gene_cluster_id) %>%
    group_by(genome_id, contig_id, contig_length) %>%
    count(name = "n_genes") %>%
    mutate(contig_length = contig_length/10^6, n_genes = n_genes/10^3) %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
    geom_point(aes(x = contig_length, y = n_genes), shape = 21, size = 2) +
    scale_x_continuous(limits = c(0, 5)) +
    scale_y_continuous(limits = c(0, 5)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "contig length (Mbp)", y = "# of gene clusters (k)")
ggsave(paste0(folder_data, "temp/17b-01-contig_length_vs_genes.png"), p, width = 4, height = 4)


# 2. plot the trees by gene content ----
# Genome tree
egcalls_wide <- egcalls %>%
    distinct(genome_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    pivot_wider(id_cols = genome_id, names_from = gene_cluster_id, values_from = value, values_fill = 0)
isolates_label <- egcalls_wide[,1] %>%
    left_join(isolates_mash) %>%
    mutate(id = genome_id) %>%
    select(id, everything())

egcc_m <- as.matrix(egcalls_wide[,-1])
dim(egcc_m)
rownames(egcc_m) <- as.character(isolates_label$genome_id)
jdm <- proxy::dist(egcc_m, method = "Jaccard")
te <- as.phylo(hclust(jdm))

p1 <- te %>%
    ggtree() %<+% isolates_label +
    geom_tippoint(aes(color = rhizobia_population), size = 2) +
    geom_tiplab(aes(label = genome_id, color = rhizobia_population), offset = 0.01) +
    geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) +
    scale_color_manual(values = rhizobia_population_colors) +
    scale_x_continuous(limits = c(-0.1, 0.6)) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "genome gene content")
p2 <- te %>%
    ggtree() %<+% isolates_label +
    geom_tippoint(aes(color = species_name), size = 2) +
    geom_tiplab(aes(label = genome_id, color = species_name), offset = 0.01) +
    geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) +
    scale_color_manual(values = ensifer_sp_colors) +
    scale_x_continuous(limits = c(-0.1, 0.6)) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "genome gene content")
p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/17b-02-tree_genomes.png"), p, width = 12, height = 6)



# 3. plot the contig tree by site and species----
egcalls_wide <- egcalls %>%
    distinct(genome_id, contig_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    pivot_wider(id_cols = c(genome_id, contig_id), names_from = gene_cluster_id, values_from = value, values_fill = 0) %>%
    unite(col = "contig_unique_id", genome_id, contig_id, remove = F)

contigs_label <- egcalls_wide[,1] %>%
    left_join(contigs_mash) %>%
    mutate(id = contig_unique_id) %>%
    select(id, everything()) %>%
    left_join(isolates_label[,-1], by = "genome_id")

egcc_m <- as.matrix(egcalls_wide[,-c(1,2,3)])
dim(egcc_m)
rownames(egcc_m) <- contigs_label$id
jdm <- proxy::dist(egcc_m, method = "Jaccard")
te <- as.dendrogram(hclust(jdm))

p1 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = rhizobia_population), size = 2) +
    geom_tiplab(aes(label = genome_id, color = rhizobia_population), offset = 0.01) +
    scale_color_manual(values = rhizobia_population_colors) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    theme_tree2(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "contig gene content")
p2 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = species_name), size = 2) +
    geom_tiplab(aes(label = genome_id, color = species_name), offset = 0.01) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    scale_color_npg() +
    theme_tree2(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "contig gene content")

p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/17b-03-tree_contigs.png"), p, width = 15, height = 15)


# 4. plot the contig tree contig type----
p1 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = ge_name), size = 2) +
    geom_tiplab(aes(label = ge_name, color = ge_name), offset = 0.01) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    theme_tree2(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""), ncol = 2)) +
    labs(title = "contig gene content")
p2 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = ge_type), size = 2) +
    geom_tiplab(aes(label = ge_name, color = ge_type), offset = 0.01) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    scale_color_npg() +
    theme_tree2(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "contig gene content")

p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/17b-04-tree_contigs.png"), p, width = 15, height = 15)














