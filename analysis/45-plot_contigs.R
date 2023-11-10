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

# 0.5 read the raw data



# 1. plot the gene numbers on each contig
p <- egcc %>%
    distinct(gene_cluster_id, gene_callers_id, .keep_all = T) %>%
    mutate(contig = str_remove_all(contig, "0")) %>%
    group_by(genome_name, contig) %>%
    summarize(n_genes = n()) %>%
    # Remove very small contigs
    filter(n_genes > 10) %>%
    # Rename the contigs by size
    arrange(genome_name, desc(n_genes)) %>%
    mutate(contig = paste0("contig_", 1:n())) %>%
    mutate(n_genes = n_genes / 1000) %>%
    ggplot() +
    geom_col(aes(x = genome_name, y = n_genes, fill = contig), color = "black", position = position_stack(reverse = T)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 9, "Set1")) +
    scale_y_continuous(breaks = seq(0, 8, 2), minor_breaks = 1:10) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(linetype = 1, linewidth = 1),
        panel.grid.minor.y = element_line(linetype = 2, linewidth = 0.3)
    ) +
    guides(fill = guide_legend(ncol = 1, reverse = T)) +
    labs(x = "genome", y = "# of genes (k)")

ggsave(paste0(folder_data, "temp/45-01-contig_genes.png"), p, width = 6, height = 4)























