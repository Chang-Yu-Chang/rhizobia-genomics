#' This script runs population genetics analysis

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
library(ggtree)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F)
isolates_mash <- read_csv(paste0(folder_data, "temp/14-isolates_mash.csv"), show_col_types = F)
genome_kmer <- read_delim(paste0(folder_data, "genomics/popgen/genome_kmer/genome_kmer.txt"), show_col_types = F)
list_sig <- read_delim(paste0(folder_data, "genomics/popgen/genome_kmer/list_sig.txt"), delim = "\t", col_names = "file_name")
# isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"), show_col_types = F)
# genomes_mapping <- read_csv(paste0(folder_data, "temp/00-genomes_mapping.csv"), show_col_types = F)


# 0. Clean uo ---
# 0.1 SNPs ----
snp <- read_table(paste0(folder_data, "genomics/popgen/usda1106/snippy/core.tab"), show_col_types = F)
gs <- read_table(paste0(folder_data, "genomics/popgen/usda1106/snippy/core.txt"), show_col_types = F)
vcf <- read.vcfR(paste0(folder_data, "genomics/popgen/usda1106/snippy/core.vcf"))

# SNPs on chromosome
vcfR::getCHROM(vcf) %>% table()
# VCF to genlight
gl <- vcfR2genlight(vcf)
# Rearrange the metadata
isolates_arr <- isolates %>%
    filter(!is.na(exp_id), exp_id != "ncbi") %>%
    left_join(isolates_mash) %>%
    arrange(match(genome_id, gl$ind.names))
# PCA
pca <- glPca(gl, nf = 10, center = T, scale = T)
# convert scores of vcf.pca into a tibble
sc <- as_tibble(pca$scores) %>%
    bind_cols(isolates_arr)

# 0.2 k-mer ----
list_sig <- list_sig %>%
    mutate(genome_name = str_remove(file_name, ".*/genomes/") %>% str_remove("/04-taxonomy.*")) %>%
    left_join(isolates_mash)

genome_kmer <- genome_kmer %>%
    mutate(row_name = colnames(.)) %>%
    pivot_longer(-row_name, names_to = "col_name", values_to = "distance") %>%
    # Clean the matrix name
    mutate(row_name = str_remove(row_name, ".*/genomes/") %>% str_remove("/02.*")) %>%
    mutate(col_name = str_remove(col_name, ".*/genomes/") %>% str_remove("/02.*")) %>%
    pivot_wider(names_from = col_name, values_from = distance)

genome_kmer_long <- genome_kmer %>%
    pivot_longer(-row_name, names_to = "col_name", values_to = "distance") %>%
    left_join(rename_with(list_sig, function (x) paste0(x, 1)) %>% rename(col_name = genome_name1)) %>%
    left_join(rename_with(list_sig, function (x) paste0(x, 2)) %>% rename(row_name = genome_name2))


# 1. PCs ----
pcs <- round(100 * pca$eig / sum(pca$eig), 2)
p <- tibble(eigenvalue = 1:length(pca$eig), pvar = pcs) %>%
    ggplot() +
    geom_col(aes(x = eigenvalue, y = pvar)) +
    scale_x_continuous(limits = c(0.5,10), breaks = 1:10) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "eigenvalues", y = "variance explained (%)")

ggsave(paste0(folder_data, "temp/18-01-snp_PCs.png"), p, width = 4, height = 3)

# 2. PCA based on SNPs ----
set.seed(1)
p1 <- sc %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = rhizobia_population), shape = 21, size = 3, stroke = 1, position = position_jitter(width = 10, height = 10)) +
    scale_color_manual(values = rhizobia_population_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 variance = ", pcs[1], "%"), y = paste0("PC2 variance = ", pcs[2], "%"))

p2 <- sc %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = species_name), shape = 21, size = 3, stroke = 1, position = position_jitter(width = 10, height = 10)) +
    scale_color_npg() +
    scale_color_manual(values = ensifer_sp_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 variance = ", pcs[1], "%"), y = paste0("PC2 variance = ", pcs[2], "%"))

p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"))
ggsave(paste0(folder_data, "temp/18-02-snp_PCA.png"), p, width = 10, height = 3)

# 3. PCA based on medicae and meliloti only ----
ids_kept <- which(isolates_arr$species_name %in% c("Ensifer meliloti", "Ensifer medicae"))
gl_mm <- gl[ids_kept]
pca <- glPca(gl_mm, nf = 10, center = T, scale = T)
# convert scores of vcf.pca into a tibble
sc <- as_tibble(pca$scores) %>%
    bind_cols(isolates_arr[ids_kept,])
pcs <- round(100 * pca$eig / sum(pca$eig), 2)

set.seed(1)
p1 <- sc %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = rhizobia_population), shape = 21, size = 3, stroke = 1) +
    scale_color_manual(values = rhizobia_population_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 variance = ", pcs[1], "%"), y = paste0("PC2 variance = ", pcs[2], "%"))

p2 <- sc %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = species_name), shape = 21, size = 3, stroke = 1) +
    scale_color_manual(values = ensifer_sp_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 variance = ", pcs[1], "%"), y = paste0("PC2 variance = ", pcs[2], "%"))

p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"))
ggsave(paste0(folder_data, "temp/18-03-snp_PCA_mm.png"), p, width = 10, height = 3)



# 4. tree based on all 32 genomes ----
# Generated pairwise distances between samples that we will plot in a tree format
tree_data <- aboot(gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50)

tree_tips <- isolates_arr %>%
    mutate(id = genome_id) %>%
    select(id, everything())

p1 <- tree_data %>%
    ggtree() %<+% tree_tips +
    geom_tippoint(aes(color = rhizobia_population), size = 2) +
    geom_tiplab(aes(label = genome_id, color = rhizobia_population), offset = 0.01) +
    scale_color_manual(values = rhizobia_population_colors) +
    scale_x_continuous(limits = c(-0.01, 0.32)) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "SNPs")

p2 <- tree_data %>%
    ggtree() %<+% tree_tips +
    geom_tippoint(aes(color = species_name), size = 2) +
    geom_tiplab(aes(label = genome_id, color = species_name), offset = 0.01) +
    scale_color_manual(values = ensifer_sp_colors) +
    scale_x_continuous(limits = c(-0.01, 0.32)) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "SNPs")


p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/18-04-snp_tree.png"), p, width = 10, height = 6)


# 5. tree based on medicae and meliloti only ----
tree_data <- aboot(gl_mm, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50)

tree_tips <- isolates_arr[ids_kept,] %>%
    mutate(id = genome_id) %>%
    select(id, everything())

p1 <- tree_data %>%
    ggtree() %<+% tree_tips +
    geom_tippoint(aes(color = rhizobia_population), size = 2) +
    geom_tiplab(aes(label = genome_id, color = rhizobia_population), offset = 0.0001) +
    #geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) +
    scale_color_manual(values = rhizobia_population_colors) +
    scale_x_continuous(limits = c(-0.001, 0.006)) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "SNPs")

p2 <- tree_data %>%
    ggtree() %<+% tree_tips +
    geom_tippoint(aes(color = species_name), size = 2) +
    geom_tiplab(aes(label = genome_id, color = species_name), offset = 0.0001) +
    #geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) +
    scale_color_manual(values = ensifer_sp_colors) +
    scale_x_continuous(limits = c(-0.001, 0.006)) +
    theme_tree(legend.position = 'centre') +
    theme(
        legend.position = c(0.2,0.8)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "SNPs")


p <- plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A", "B"), scale = 0.9) + paint_white_background()
ggsave(paste0(folder_data, "temp/18-05-snp_tree_mm.png"), p, width = 10, height = 6)


#stat_ellipse(level = 0.95, size = 1) +


if (FALSE) {
    vcf_usda <- readVcf(paste0(folder_data, "genomics/popgen/snippy_usda1106/snippy_usda1106.vcf"))
vcf_wsm <- read.vcfR(paste0(folder_data, "genomics/popgen/snippy_wsm419/snippy_wsm419.vcf"))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F)
isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"), show_col_types = F)

# 1. VCF description ----
# number of snps on each contig
vcfR::getCHROM(vcf_usda) %>% table()
vcfR::getCHROM(vcf_wsm) %>% table()

# extract only the chromosome

as_tibble(vcf_usda@fix) %>%
    distinct(POS)
snp_usda <- as_tibble(vcf_usda@gt, .name_repair = "minimal") %>%
    setNames(c("snp", paste0("g", 1:32)))



snp_usda %>%
    clean_names() %>%
    slice(1:10000) %>%
    pivot_longer(cols = -snp) %>%
    filter(value != 1)


# 1. descriptions ----

# Convert the vcf to a genclone object
gi_usda <- vcfR2genind(snp_usda)
class(gi_usda)
gc_usda <- as.genclone(gi_usda)
class(gc_usda)

gc_usda$pop <- factor(isolates_mapping$rhizobia_population)


poppr(gc_usda)


# mll(geni_usda)
#
# mll(monpop) <- "original"
# diversity_stats(mll(monpop))
#mlg.filter(x, distance = xdis) <- 1 + .Machine$double.eps^0.5



# Remove uninformative loci
#geni_usda_cut <- informloci(geni_usda, cutoff = 2/nInd(geni_usda), MAF = 0.01, quiet = FALSE)
popNames(core_genind)


genlight_usda <- dartR::gi2gl(geni_usda)












pca_usda <- glPca(genlight_usda, nf = 2, center = T, scale = T)
eig_usda <- sort(pca_usda$eig/sum(pca_usda$eig), decreasing = T)

# Merge the tb
isolates_pca <- pca_usda$scores %>%
    as_tibble() %>%
    mutate(PC1 = PC1/100, PC2 = PC2/100) %>%
    mutate(genome_name = genlight_usda$ind.names) %>%
    left_join(isolates_ensifer)

set.seed(1)
p1 <- isolates_pca %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site), shape = 21, size = 3, stroke = 1.5, position = position_jitter(width = 0.02, height = 0.02)) +
    scale_color_manual(values = rhizobia_population_colors) +
    theme_classic() +
    theme(
        #legend.background = element_rect(fill = "white", color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = 2),
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(eig[1]*100,1), "%)"),
         y = paste0("PC2 (", round(eig[2]*100,1), "%)"))
}
















# 6. tree based on kmer matrix ----
library(ape)
m <- 1-as.matrix(genome_kmer[,-1]) # The genome_kmer is similarity -> convert it to dissimilarity
row_sums <- rowSums(m)
col_sums <- colSums(m)
m_ordered <- m[order(row_sums, decreasing = TRUE), order(col_sums, decreasing = TRUE)]
tree <- hclust(as.dist(m_ordered), method = "ward.D") %>% as.phylo()
tree_tips <- tibble(genome_name = genome_kmer$row_name) %>%
    left_join(isolates) %>%
    left_join(isolates_mash) %>%
    mutate(id = genome_name) %>%
    select(id, everything()) %>%
    mutate(genome_name = factor(genome_name, tree$tip.label)) %>%
    arrange(genome_name)

p <- tree %>%
    ggtree(ladderize = F) %<+% tree_tips +
    geom_tippoint(aes(color = species_name), size = 2) +
    geom_tiplab(aes(label = genome_id, color = species_name), offset = 0.1) +
    #geom_nodelab(size = 2, nudge_x = -0.006, nudge_y = 1) +
    scale_color_manual(values = ensifer_sp_colors) +
    scale_x_continuous(limits = c(-0.001, 5.5)) +
    theme_bw() +
    #theme_tree(legend.position = 'centre') +
    theme(
#        legend.position = c(0.2,0.8)
    ) +
    guides(color = guide_legend(override.aes = aes(label = ""))) +
    labs(title = "genome k-mer")


ggsave(paste0(folder_data, "temp/18-06-genome_kmer_tree.png"), p, width = 6, height = 5)

# 7. plot the kmer matrix ----
p <- genome_kmer_long %>%
    mutate(row_name = factor(row_name, rev(tree$tip.label)), col_name = factor(col_name, tree$tip.label)) %>%
    mutate(genome_id1 = factor(genome_id1, tree_tips$genome_id), genome_id2 = factor(genome_id2, rev(tree_tips$genome_id))) %>%
    ggplot() +
    geom_tile(aes(x = genome_id1, y = genome_id2, fill = distance)) +
    #geom_tile(aes(x = row_name, y = col_name, fill = distance)) +
    scale_fill_gradient(low = "snow", high = "maroon") +
    scale_x_discrete(position = "top") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = -1),
        axis.title = element_blank()
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/18-07-genome_kmer.png"), p, width = 6, height = 5)








