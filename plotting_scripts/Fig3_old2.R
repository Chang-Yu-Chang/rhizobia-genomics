#' This script plots the network based one distances

renv::load()
library(tidyverse)
library(cowplot)
# library(tidygraph)
# library(ggraph)
# library(ggforce)
library(ape)
library(tidytree)
library(ggtree)
library(proxy) # For computing jaccard distance
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
library(vegan) # for computing jaccard
source(here::here("metadata.R"))

gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
species_shapes <- c(meliloti = 21, medicae = 22, adhaerens = 15, canadensis = 16)
load(file = paste0(folder_data, "phylogenomics_analysis/networks/networks.rdata"))
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

if (FALSE) {

# Panel A. kmer networks
thr_kmer <- 0.85
p1 <- g %>%
    activate(edges) %>%
    mutate(kmer_bin = ifelse(d_kmer < thr_kmer, paste0("d_kmer<", thr_kmer), "nope")) %>%
    ggraph(layout = "linear", circular = T) +
    #geom_mark_hull(aes(x, y, group = species, fill = species), concavity = 4, expand = unit(3, "mm"), alpha = 0.25) +
    geom_edge_arc(aes(colour = kmer_bin), alpha = 0.1) +
    geom_node_point(aes(fill = species), shape = 21, color = "black", size = 3) +
    scale_edge_color_manual(values = c("black", "NA")) +
    scale_fill_manual(values = species_colors) +
    theme_void() +
    theme(
        plot.background = element_rect(color = NA, fill = "white"),
        plot.margin = unit(c(1,1,1,1), "mm"),
        legend.position = "right"
    ) +
    guides(edge_colour = "none") +
    labs()
}

# Panel A. PAP heatmap
p1 <- gpatl %>%
    ggplot() +
    geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 5, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome")


# Panel B. matched tree
# Plot core tree
p2_1 <- tr %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    labs()

# Plot accessory tree
p2_2 <- tr_acce %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label, color = species), hjust = 0) +
    scale_color_manual(values = species_colors)

d1 <- p2_1$data
d2 <- p2_2$data
## reverse x-axis and set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
dd <- bind_rows(d1, d2) %>% filter(isTip)

p2 <- p2_1 + geom_tree(data = d2) +
    ggnewscale::new_scale_fill() +
    geom_line(data = dd, aes(x, y, group = label), color = "grey90", linetype = 1, linewidth = .2) +
    geom_tippoint(aes(color = species), size = 2) +
    geom_tippoint(data = d2, aes(color = species), size = 2) +
    scale_x_continuous(limits = c(-0.5, 3)) +
    annotate("text", x = c(-0.3, 2.5), y = c(30,30), label = c("core genes", "gene content")) +
    theme_tree() +
    theme(
        legend.position = "top"
    )
    #theme_classic() +
    labs()



# Panel C. SNPs
snp <- read_table(paste0(folder_data, "genomics/variants/em1021/core.tab"))
gs <- read_table(paste0(folder_data, "genomics/variants/em1021/core.txt"))
vcf <- read.vcfR(paste0(folder_data, "genomics/variants/em1021/core.vcf"))

vcfR::getCHROM(vcf) %>% table() # SNPs on chromosome
gl <- vcfR2genlight(vcf) # VCF to genlight
snps <- tab(gl); dim(snps) # Number of SNPs 18664

## Test VA populations: high vs low elevation
isolates_i1 <- isolates %>% filter(population == "VA")
m <- snps[match(isolates_i1$genome_id, gl@ind.names),]; dim(m) # Filter the snps table for focal genomes
#table(apply(m, 2, sum))
dist_m1 <- vegdist(m, "jaccard")
permanova1 <- adonis2(dist_m1 ~ site_group, data = isolates_i1)
permanova1 # no
pcoa1 <- cmdscale(dist_m1, eig = T)
eigs1 <- round(pcoa1$eig / sum(pcoa1$eig)*100, 2)

## Test PA populations: urban vs suburban elevation
isolates_i2 <- isolates %>% filter(population == "PA")
m <- snps[match(isolates_i2$genome_id, gl@ind.names),]; dim(m) # Filter the snps table for focal genomes
#m <- m[,apply(m, 2, sum) != 0]; dim(m) # 323 snps??
dist_m2 <- vegdist(m, "jaccard")
permanova2 <- adonis2(dist_m2 ~ site_group, data = isolates_i2)
permanova2 # no
pcoa2 <- cmdscale(dist_m2, eig = T)
eigs2 <- round(pcoa2$eig / sum(pcoa2$eig)*100, 2)

## Plot the mds
p3_1 <- isolates_i1 %>%
    left_join(isolates_contigs) %>%
    bind_cols(tibble(mds1 = pcoa1$points[,1], mds2 = pcoa1$points[,2])) %>%
    drop_na(species) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group, shape = species), size = 3, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    scale_shape_manual(values = species_shapes) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs1[1], "%)"), y = paste0("PCoA axis 1(", eigs1[2], "%)"))

p3_2 <- isolates_i2 %>%
    left_join(isolates_contigs) %>%
    bind_cols(tibble(mds1 = pcoa2$points[,1], mds2 = pcoa2$points[,2])) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group, shape = species), size = 3, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    scale_shape_manual(values = species_shapes) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",

        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs2[1], "%)"), y = paste0("PCoA axis 1(", eigs2[2], "%)"))


# Panel D. GCVs
isolates_test <- isolates %>% filter(genome_id %in% gpat$genome_id, !is.na(site_group))
gpa_test <- gpat %>% filter(genome_id %in% isolates_test$genome_id)
dim(gpa_test) # 26887 gene clusters

## Test VA populations: high vs low elevation
isolates_i1 <- isolates_test %>% filter(population == "VA")
labels <- isolates_i1$site_group; labels <- as.factor(labels)
gpa_i <- gpa_test %>% filter(genome_id %in% isolates_i1$genome_id)
m <- as.matrix(gpa_i[,-1]); dim(m)
m <- m[,apply(m, 2, sum) != 0]; dim(m) # 23575 genes
dist_m1 <- vegdist(m, "jaccard")
permanova1 <- adonis2(dist_m1 ~ site_group, data = isolates_i1)
permanova1 # no
pcoa1 <- cmdscale(dist_m1, eig = T)
eigs1 <- round(pcoa1$eig / sum(pcoa1$eig)*100, 2)

## Test PA populations: urban vs suburban elevation
isolates_i2 <- isolates_test %>% filter(population == "PA")
labels <- isolates_i2$site_group; labels <- as.factor(labels)
gpa_i <- gpa_test %>% filter(genome_id %in% isolates_i2$genome_id)
m <- as.matrix(gpa_i[,-1]); dim(m)
m <- m[,apply(m, 2, sum) != 0]; dim(m) # 12946 genes
dist_m2 <- vegdist(m, "jaccard")
permanova2 <- adonis2(dist_m2 ~ site_group, data = isolates_i2)
permanova2 # no
pcoa2 <- cmdscale(dist_m2, eig = T)
eigs2 <- round(pcoa2$eig / sum(pcoa2$eig)*100, 2)

## Plot the mds
p4_1 <- isolates_i1 %>%
    left_join(isolates_contigs) %>%
    bind_cols(tibble(mds1 = pcoa1$points[,1], mds2 = pcoa1$points[,2])) %>%
    drop_na(species) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group, shape = species), size = 3, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    scale_shape_manual(values = species_shapes) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",
        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs1[1], "%)"), y = paste0("PCoA axis 1(", eigs1[2], "%)"))

p4_2 <- isolates_i2 %>%
    left_join(isolates_contigs) %>%
    bind_cols(tibble(mds1 = pcoa2$points[,1], mds2 = pcoa2$points[,2])) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group, shape = species), size = 3, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",
        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs2[1], "%)"), y = paste0("PCoA axis 1(", eigs2[2], "%)"))


p_left <- plot_grid(p1 + ggtitle("Gene presence/absence"), p2, nrow = 2, scale = 0.9, labels = c("A", "B"))
leg1 <- cowplot::get_plot_component(p3_1 + theme(legend.background = element_rect(color = "black")), 'guide-box-top', return_all = TRUE)
leg2 <- cowplot::get_plot_component(p3_2 + theme(legend.background = element_rect(color = "black")), 'guide-box-top', return_all = TRUE)
p_right <- plot_grid(
    p3_1 + guides(color = "none"), p3_2 + guides(color = "none"),
    p4_1, p4_2, nrow = 2, rel_heights = c(1, 1), scale = 0.95,
    align = "v", axis = "lr", labels = c("C", "", "D", ""))


#p_top <- plot_grid(p1 + ggtitle("Jaccard distance"), p2 + ggtitle("Gene presence/absence"), nrow = 1, labels = c("A", "B"), scale = 0.8)
p <- plot_grid(p_left, p_right, nrow = 1, labels = c("", ""), scale = c(1, 0.9), rel_heights = c(1, 1.2)) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 6)




