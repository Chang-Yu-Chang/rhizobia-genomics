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


# 2. plot the gene content tree ----
egcalls_wide <- egcalls %>%
    distinct(genome_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    pivot_wider(id_cols = genome_id, names_from = gene_cluster_id, values_from = value, values_fill = 0)
isolates_label <- egcalls_wide[,1] %>% left_join(isolates_mash)

egcc_m <- as.matrix(egcalls_wide[,-1])
dim(egcc_m)
rownames(egcc_m) <- as.character(isolates_label$genome_id)
jdm <- proxy::dist(egcc_m, method = "Jaccard")
te <- as.phylo(hclust(jdm))
#clade_highlight <- tibble(node = c(109, 114, 116), ge_type = c("chromosome", "pSymA like", "pSymB like"))
te$tip.label
isolates_label <- tibble(node = 1:length(te$tip.label), genome_id = te$tip.label) %>% left_join(isolates_label)

table_data <- data.frame(tip_label = c(1, 2, 3, 4),
                         label = c("Label_A", "Label_B", "Label_C", "Label_D"))

merged_data <- left_join(data.frame(tip_label = tip_labels), table_data, by = c("tip_label"))

p <- te %>%
    #full_join(rename(isolates_label, tip.label = genome_id), by = "tip.label") %>%
    ggtree(layout = "rectangular") +
    #layout_circular() +
    geom_tiplab() +
    geom_tiplab(data = isolates_label, aes(label = genome_id)) +
    #geom_text(aes(label=node), hjust=-.3) +
    # geom_hilight(data = clade_highlight, aes(node = node, fill = ge_type), type = "roundrect") +
    theme_tree() +
    theme(
        plot.margin = unit(c(10, 10, 10, 10), "mm")
    ) +
    labs()

ggsave(paste0(folder_data, "temp/17b-02-tree_genomes.png"), p, width = 10, height = 10)



# 3. plot the contig tree ----
egcalls_wide <- egcalls %>%
    distinct(genome_id, contig_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    pivot_wider(id_cols = c(genome_id, contig_id), names_from = gene_cluster_id, values_from = value, values_fill = 0) %>%
    unite(col = "contig_unique_id", genome_id, contig_id, remove = F)

contigs_label <- egcalls_wide[,1] %>% left_join(contigs_mash)
egcc_m <- as.matrix(egcalls_wide[,-c(1,2,3)])
dim(egcc_m)
rownames(egcc_m) <- as.character(contigs_label$ge_name)
jdm <- proxy::dist(egcc_m, method = "Jaccard")
te <- as.dendrogram(hclust(jdm))
clade_highlight <- tibble(node = c(109, 114, 116), ge_type = c("chromosome", "pSymA like", "pSymB like"))

p <- te %>%
    ggtree(layout = "rectangular") +
    layout_circular() +
    geom_tiplab() +
    #geom_text(aes(label=node), hjust=-.3) +
    geom_hilight(data = clade_highlight, aes(node = node, fill = ge_type), type = "roundrect") +
    theme(
        plot.margin = unit(c(10, 10, 10, 10), "mm")
    ) +
    labs()

ggsave(paste0(folder_data, "temp/17b-03-tree_contigs.png"), p, width = 15, height = 15)

















