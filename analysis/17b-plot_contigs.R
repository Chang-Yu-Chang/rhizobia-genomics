#' This script plots the gene content by the contigs

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
library(proxy) # for computing jaccard matrix
library(ape)
library(ggtree)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    filter(!genome_id %in% paste0("g", c(1,7,12,14,18)))
egc <- read_csv(paste0(folder_data, "temp/17-egc.csv"), show_col_types = F)
gcalls <- read_csv(paste0(folder_data, "temp/17-gcalls.csv"), show_col_types = F)
egcalls <- read_csv(paste0(folder_data, "temp/17-egcalls.csv"), show_col_types = F)
nrow(isolates) # This should be 32 ensifer isolates + 4 ncbi genomes

# 0. Data wrangling ----


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



# 3. plot the contig tree by site and mash species----
egcalls_wide <- egcalls %>%
    distinct(genome_id, contig_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    pivot_wider(id_cols = c(genome_id, contig_id), names_from = gene_cluster_id, values_from = value, values_fill = 0) %>%
    unite(col = "contig_unique_id", genome_id, contig_id, remove = F)

contigs_label <- egcalls_wide[,1] %>%
    left_join(contigs_mash) %>%
    left_join(contigs_blast) %>%
    mutate(id = contig_unique_id) %>%
    select(id, everything()) %>%
    left_join(isolates_label[,c("genome_id", "rhizobia_population")], by = "genome_id")

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
    theme_tree(legend.position = 'centre') +
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
    scale_color_manual(values = ensifer_sp_colors) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "contig gene content; blast species")

p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/17b-03-tree_contigs.png"), p, width = 15, height = 15)


# 4. plot the contig tree by mash contig type----
p1 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = ge_name), size = 2) +
    geom_tiplab(aes(label = ge_name, color = ge_name), offset = 0.01) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""), ncol = 2)) +
    labs(title = "contig gene content; mash")
p2 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = ge_type), size = 2) +
    geom_tiplab(aes(label = ge_name, color = ge_type), offset = 0.01) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    scale_color_npg() +
    #scale_color_manual(values = ensifer_sp_colors) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "contig gene content; mash")

p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/17b-04-tree_contigs_mash.png"), p, width = 15, height = 15)


# 5. plot the contig tree by blast contig type----
p1 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = contig_name), size = 2) +
    geom_tiplab(aes(label = contig_name, color = contig_name), offset = 0.01) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""), ncol = 2)) +
    labs(title = "contig gene content; blast")
p2 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = species_name), size = 2) +
    geom_tiplab(aes(label = contig_name, color = species_name), offset = 0.01) +
    scale_x_continuous(limits = c(-1.1, 0.3)) +
    scale_color_manual(values = ensifer_sp_colors) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.9)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "contig gene content; blast")

p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/17b-05-tree_contigs_blast.png"), p, width = 15, height = 15)




















# 6. plot the genome tree with heatmap ----
#library(aplot) # For aligning tree to the heatmap
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

# Tree
p1 <- te %>%
    ggtree() %<+% isolates_label +
    geom_tippoint(aes(color = species_name), size = 3) +
    geom_tiplab(aes(label = genome_id, fill = species_name), offset = 0.01) +
    #geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) +
    #scale_color_manual(values = rhizobia_population_colors) +
    scale_color_manual(values = ensifer_sp_colors[1:2]) +
    scale_x_continuous(limits = c(-0.05, 0.55)) +
    scale_y_continuous(limits = c(0.5, 36.5), expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    #theme_classic() +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        #plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""), title = "species by mash")) +
    labs(title = "")

# heatmap of species number
p2 <- isolates_label %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    arrange(genome_id) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = 1, fill = rhizobia_population), width = 0.9) +
    scale_fill_manual(values = rhizobia_population_colors) +
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    theme_tree() +
    theme(
        legend.position = "top",
        legend.key.height = unit(1, "mm"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
    ) +
    guides(fill = guide_legend(ncol = 1, title = NULL, byrow = T)) +
    labs(x = NULL, y = NULL)

# heatmap of species name
p3 <- isolates_label %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    arrange(genome_id) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = 1, fill = species_name), width = 0.9) +
    scale_fill_manual(values = ensifer_sp_colors[1:2], labels = c("medicae", "meliloti")) +
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    theme_tree() +
    theme(
        legend.position = "top",
        legend.key.height = unit(1, "mm"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
    ) +
    guides(fill = guide_legend(ncol = 1, title = NULL, byrow = T)) +
    labs(x = NULL, y = NULL)

# heatmap of growth trait
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'), show_col_types = F)
isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    group_by(exp_id) %>%
    arrange(exp_id, temperature) %>%
    #pivot_wider(names_from = temperature, values_from = c(r, lag, maxOD)) %>%
    left_join(isolates_label) %>%
    ungroup() %>%
    bind_rows(
        tibble(genome_id = rep(c("em1021", "em1022", "usda1106", "wsm419"), each = 4),
               temperature = rep(c("25c", "30c", "35c", "40c"), 4),
               r = rep(NA, 16),
               lag = rep(NA, 16),
               maxOD = rep(NA, 16))
    )

median_r <- median(isolates_gc$r, na.rm = T)
median_lag <- median(isolates_gc$lag, na.rm = T)
median_maxOD <- median(isolates_gc$maxOD, na.rm = T)


# Growth rate
p4 <- isolates_gc  %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    arrange(genome_id) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = temperature, fill = r), width = 0.9) +
    #scale_fill_manual(values = ensifer_sp_colors[1:2]) +
    scale_fill_gradient2(high = "maroon", low = "steelblue", mid = "snow", midpoint = median_r, breaks = rev(c(0.1, 0.5, 1, 1.5)), limits = c(0,1.5)) +
    #scale_fill_gradient2() +
    scale_x_discrete(expand = c(0,0), drop = F) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    theme_minimal() +
    #theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
    ) +
    guides(fill = guide_colourbar(title = "r")) +
    labs(x = NULL, y = NULL)

# Lag time
p5 <- isolates_gc  %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    arrange(genome_id) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = temperature, fill = lag), width = 0.9) +
    scale_fill_gradient2(high = "steelblue", low = "maroon", mid = "snow", midpoint = median_lag, breaks = rev(seq(12,48,12)), limits = c(0,40)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    theme_minimal() +
    #theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
    ) +
    guides(fill = guide_colorbar(title = "lag")) +
    labs(x = NULL, y = NULL)

# maxOD
p6 <- isolates_gc  %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    arrange(genome_id) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = temperature, fill = maxOD), width = 0.9) +
    scale_fill_gradient2(high = "maroon", low = "steelblue", mid = "snow", midpoint = median_maxOD, breaks = rev(seq(0.1, 0.5, 0.1))) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    theme_minimal() +
    #theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
    ) +
    guides(fill = guide_colorbar(title = "maxOD")) +
    labs(x = NULL, y = NULL)

p <- plot_grid(p1,p2,p4,p5,p6, rel_widths = c(1,0.15,0.5,0.5,0.5),
               nrow = 1, align = "h", axis = "tb") + paint_white_background()

ggsave(paste0(folder_data, "temp/17b-06-tree_genome_heatmap.png"), p, width = 16, height = 10)

# 7. plot the contig tree with heatmap ----
egcalls_wide <- egcalls %>%
    distinct(genome_id, contig_id, gene_cluster_id) %>%
    mutate(value = 1) %>%
    pivot_wider(id_cols = c(genome_id, contig_id), names_from = gene_cluster_id, values_from = value, values_fill = 0) %>%
    unite(col = "contig_unique_id", genome_id, contig_id, remove = F)

contigs_label <- egcalls_wide[,1] %>%
    left_join(contigs_mash) %>%
    left_join(contigs_blast) %>%
    # Clean the unique contig id
    mutate(contig_unique_id = str_remove(contig_unique_id, "0000000000")) %>%
    mutate(id = contig_unique_id) %>%
    select(id, everything()) %>%
    left_join(isolates_label[,c("genome_id", "rhizobia_population")], by = "genome_id")

egcc_m <- as.matrix(egcalls_wide[,-c(1,2,3)])
dim(egcc_m)
rownames(egcc_m) <- contigs_label$id
jdm <- proxy::dist(egcc_m, method = "Jaccard")
te <- as.dendrogram(hclust(jdm))

# egcalls_wide <- egcalls %>%
#     distinct(genome_id, gene_cluster_id) %>%
#     mutate(value = 1) %>%
#     pivot_wider(id_cols = genome_id, names_from = gene_cluster_id, values_from = value, values_fill = 0)
# isolates_label <- egcalls_wide[,1] %>%
#     left_join(isolates_mash) %>%
#     mutate(id = genome_id) %>%
#     select(id, everything())
#
# egcc_m <- as.matrix(egcalls_wide[,-1])
# dim(egcc_m)
# rownames(egcc_m) <- as.character(isolates_label$genome_id)
# jdm <- proxy::dist(egcc_m, method = "Jaccard")
# te <- as.phylo(hclust(jdm))

# Tree
p1 <- te %>%
    ggtree() %<+% contigs_label +
    geom_tippoint(aes(color = species_name), size = 3) +
    geom_tiplab(aes(label = contig_unique_id, color = species_name), offset = 0.01, align = T) +
    #geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) +
    #scale_color_manual(values = rhizobia_population_colors) +
    scale_color_manual(values = ensifer_sp_colors[1:2]) +
    # scale_x_continuous(limits = c(-0.05, 0.55)) +
    # scale_y_continuous(limits = c(0.5, 36.5), expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    #theme_classic() +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8),
        plot.margin = unit(c(0, 20, 0, 0), "mm"),
        #plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""), title = "species by mash")) +
    labs(title = "")
p1
# heatmap of spcies population
p2 <- contigs_label %>%
    mutate(contig_unique_id = factor(contig_unique_id, rev(get_taxa_name(p1)))) %>%
    arrange(contig_unique_id) %>%
    ggplot() +
    geom_tile(aes(x = contig_unique_id, y = 1, fill = rhizobia_population), width = 0.9) +
    scale_fill_manual(values = rhizobia_population_colors) +
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    theme_tree() +
    theme(
        legend.position = "top",
        legend.key.height = unit(1, "mm"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
    ) +
    guides(fill = guide_legend(ncol = 1, title = NULL, byrow = T)) +
    labs(x = NULL, y = NULL)

# heatmap of contig name by blast
p3 <- contigs_label %>%
    mutate(contig_unique_id = factor(contig_unique_id, rev(get_taxa_name(p1)))) %>%
    arrange(contig_unique_id) %>%
    ggplot() +
    geom_tile(aes(x = contig_unique_id, y = 1, fill = contig_name), width = 0.9) +
    scale_x_discrete(expand = c(0,0)) +
    scale_fill_manual(values = contig_colors) +
    coord_flip() +
    theme_tree() +
    theme(
        legend.position = "top",
        legend.key.height = unit(1, "mm"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
    ) +
    guides(fill = guide_legend(ncol = 1, title = NULL, byrow = T)) +
    labs(x = NULL, y = NULL)

# heatmap of growth trait
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'), show_col_types = F)
isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    group_by(exp_id) %>%
    arrange(exp_id, temperature) %>%
    left_join(isolates_label[,c("genome_id", "exp_id")]) %>%
    ungroup() %>%
    bind_rows(
        tibble(genome_id = rep(c("em1021", "em1022", "usda1106", "wsm419"), each = 4),
               temperature = rep(c("25c", "30c", "35c", "40c"), 4),
               r = rep(NA, 16),
               lag = rep(NA, 16),
               maxOD = rep(NA, 16))
    ) %>%
    select(-exp_id)


contigs_gc <- contigs_label %>%
    select(genome_id, contig_unique_id) %>%
    full_join(isolates_gc, relationship = "many-to-many")


# Growth rate
p4 <- contigs_gc  %>%
    mutate(contig_unique_id = factor(contig_unique_id, rev(get_taxa_name(p1)))) %>%
    arrange(contig_unique_id) %>%
    ggplot() +
    geom_tile(aes(x = contig_unique_id, y = temperature, fill = r), width = 0.9) +
    #scale_fill_manual(values = ensifer_sp_colors[1:2]) +
    scale_fill_gradient2(high = "maroon", low = "steelblue", mid = "snow", midpoint = median_r, breaks = rev(c(0.1, 0.5, 1, 1.5)), limits = c(0,1.5)) +
    #scale_fill_gradient2() +
    scale_x_discrete(expand = c(0,0), drop = F) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    theme_minimal() +
    #theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
    ) +
    guides(fill = guide_colourbar(title = "r")) +
    labs(x = NULL, y = NULL)

# Lag time
p5 <- contigs_gc  %>%
    mutate(contig_unique_id = factor(contig_unique_id, rev(get_taxa_name(p1)))) %>%
    arrange(contig_unique_id) %>%
    ggplot() +
    geom_tile(aes(x = contig_unique_id, y = temperature, fill = lag), width = 0.9) +
    scale_fill_gradient2(high = "steelblue", low = "maroon", mid = "snow", midpoint = median_lag, breaks = rev(seq(12,48,12)), limits = c(0,40)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    theme_minimal() +
    #theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
    ) +
    guides(fill = guide_colorbar(title = "lag")) +
    labs(x = NULL, y = NULL)

# maxOD
p6 <- contigs_gc  %>%
    mutate(contig_unique_id = factor(contig_unique_id, rev(get_taxa_name(p1)))) %>%
    arrange(contig_unique_id) %>%
    ggplot() +
    geom_tile(aes(x = contig_unique_id, y = temperature, fill = maxOD), width = 0.9) +
    scale_fill_gradient2(high = "maroon", low = "steelblue", mid = "snow", midpoint = median_maxOD, breaks = rev(seq(0.1, 0.5, 0.1))) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    theme_minimal() +
    #theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
    ) +
    guides(fill = guide_colorbar(title = "maxOD")) +
    labs(x = NULL, y = NULL)

p <- plot_grid(p1,p3,p2,p4,p5,p6, rel_widths = c(1,0.15,0.15,0.5,0.5,0.5),
               nrow = 1, align = "h", axis = "tb") + paint_white_background()

ggsave(paste0(folder_data, "temp/17b-07-tree_contig_heatmap.png"), p, width = 20, height = 20)


isolates_gc %>%
    left_join(isolates_label) %>%
    filter(species_name=="Ensifer medicae") %>%
    filter(temperature=="30c")












