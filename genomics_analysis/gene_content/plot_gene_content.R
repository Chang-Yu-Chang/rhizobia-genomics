#' This script plots the heatmap of gene content

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv"))
gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpatl.csv")) %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))

# 1. gene presence absence heatmap ----
p <- gpatl %>%
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

ggsave(paste0(folder_data, "genomics_analysis/gene_content/01-gpa_heatmap.png"), p, width = 6, height = 3)
nrow(gene_order) # 26504 genes in the pangenome

# 2. gene frequency spectrum ----
gfs <- gpa %>%
    pivot_longer(-gene) %>%
    group_by(gene) %>%
    summarize(n_genomes = sum(value)) %>%
    group_by(n_genomes) %>%
    count()
p <- gfs %>%
    ggplot() +
    geom_col(aes(x = n_genomes, y = n), color = "black", fill = "white") +
    scale_x_continuous(breaks = c(1, 31, seq(5, 30, 5))) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "# of genomes", y = "# of genes")

ggsave(paste0(folder_data, "genomics_analysis/gene_content/02-gene_frequency_spectrum.png"), p, width = 5, height = 4)

# 3. singletons per genome ----
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
list_sg <- gpa$gene[which(apply(gpa[,-1], 1, sum) == 1)]
gpatlsg <- gpatl %>%
    filter(gene %in% list_sg) %>%
    left_join(isolates_tax) %>%
    mutate(contig_species = str_remove(contig_species, "E. ")) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))


p1 <- gpatlsg %>%
    group_by(genome_id, contig_species) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = n, fill = contig_species), color = "black") +
    scale_fill_manual(values = species_colors) +
    theme_bw() +
    theme() +
    guides() +
    labs(title = "Symbiotic and non-symbiotic strains", y = "# of singltons")

p2 <- gpatlsg %>%
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    group_by(genome_id, contig_species) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = n, fill = contig_species), color = "black") +
    scale_fill_manual(values = species_colors) +
    theme_bw() +
    theme() +
    guides() +
    labs(title = "Symbiotic strains", y = "# of singltons")

p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
ggsave(paste0(folder_data, "genomics_analysis/gene_content/03-singletons.png"), p, width = 10, height = 8)

#
gpatlsg %>%
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    group_by(genome_id, contig_species) %>%
    count() %>%
    pull(n) %>%
    range()
















