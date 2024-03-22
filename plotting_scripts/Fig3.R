#' This script generates figure 3

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
library(vegan) # for computing jaccard
source(here::here("analysis/00-metadata.R"))

# Read gene presence-absence data
gpa <- read_csv(paste0(folder_data, "temp/13-gpa.csv"))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))

#
removed_st <- group_by(isolates_contigs, genome_id) %>% count() %>% filter(n !=1)
isolates_sp <- select(isolates_contigs, genome_id, species) %>% filter(!genome_id %in% removed_st$genome_id)
species_shapes <- c(meliloti = 21, medicae = 22, adhaerens = 15, canadensis = 16)

# Panel A. Cartoons
p1 <- ggdraw() + draw_text("placeholder")


# Panel B. pangenome
gpal <- gpa %>%
    #gpa[1:10,c(1, 3:100, 30000:31000)] %>%
    pivot_longer(-name, names_to = "gene") %>%
    mutate(gene = factor(gene, names(gpa))) %>%
    mutate(value = factor(value))

p2 <- gpal %>%
    ggplot() +
    geom_tile(aes(x = gene, y = name, fill = value)) +
    scale_fill_manual(values = c(`1` = "black", `0` = "grey95"), labels = c(`1` = "presence", `0` = "absence")) +
    coord_radial(start = 0.25 * pi, end = 1.75 * pi, inner.radius = 0.3, direction = -1) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.r = element_text(color = "red", angle = 90),
        legend.position.inside = c(0.5, 0.5)
    ) +
    guides(fill = guide_legend(position = "bottom", title = NULL)) +
    labs(y = "")


# Panel C. SNPs
snp <- read_table(paste0(folder_data, "genomics/variants/em1021_snippy/core.tab"))
gs <- read_table(paste0(folder_data, "genomics/variants/em1021_snippy/core.txt"))
vcf <- read.vcfR(paste0(folder_data, "genomics/variants/em1021_snippy/core.vcf"))

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
    left_join(isolates_sp) %>%
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
    left_join(isolates_sp) %>%
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







# Panel C. GCVs
isolates_test <- isolates %>% filter(genome_id %in% gpa$name, !is.na(site_group))
gpa_test <- gpa %>% filter(name %in% isolates_test$genome_id)
dim(gpa_test) # 31965 gene clusters

## Test VA populations: high vs low elevation
isolates_i1 <- isolates_test %>% filter(population == "VA")
labels <- isolates_i1$site_group; labels <- as.factor(labels)
gpa_i <- gpa_test %>% filter(name %in% isolates_i1$genome_id)
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
gpa_i <- gpa_test %>% filter(name %in% isolates_i2$genome_id)
m <- as.matrix(gpa_i[,-1]); dim(m)
m <- m[,apply(m, 2, sum) != 0]; dim(m) # 12946 genes
dist_m2 <- vegdist(m, "jaccard")
permanova2 <- adonis2(dist_m2 ~ site_group, data = isolates_i2)
permanova2 # no
pcoa2 <- cmdscale(dist_m2, eig = T)
eigs2 <- round(pcoa2$eig / sum(pcoa2$eig)*100, 2)

## Plot the mds
p4_1 <- isolates_i1 %>%
    left_join(isolates_sp) %>%
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
    left_join(isolates_sp) %>%
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


leg1 <- cowplot::get_plot_component(p3_1 + theme(legend.background = element_rect(color = "black")), 'guide-box-top', return_all = TRUE)
leg2 <- cowplot::get_plot_component(p3_2 + theme(legend.background = element_rect(color = "black")), 'guide-box-top', return_all = TRUE)
#p_left <- plot_grid(p1, p2, nrow = 2)
p_right <- plot_grid(
    leg1, leg2,
    p3_1 + guides(color = "none"), p3_2 + guides(color = "none"),
    p4_1, p4_2, nrow = 3, rel_heights = c(0.2, 1, 1), scale = 0.95,
    align = "v", axis = "lr", labels = c("", "", "B", "", "C", ""))

p <- plot_grid(p2, p_right, rel_widths = c(1,1), labels = c("A", "")) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 6)


# Number of genes
gpa %>%
    rename(genome_id = name) %>%
    pivot_longer(-genome_id) %>%
    group_by(name) %>%
    summarize(n_genomes = sum(value)) %>%
    filter(n_genomes < 41)
