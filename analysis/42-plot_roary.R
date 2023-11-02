#' This script plot the pangenome analysis

library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
source(here::here("analysis/00-metadata.R"))

# 0. read data
# gpa_meta <- read_csv("~/mgjw/problem_set7/roary/gene_presence_absence.csv", show_col_types = F)
# gpa <- read_table("~/mgjw/problem_set7/roary/gene_presence_absence.Rtab", show_col_types = F) %>% select(-`NC_017512.1:`)

folder_roary <- paste0(folder_data, "temp/plasmidsaurus/summary/42-roary/roary2/")
gpa_meta <- read_csv(paste0(folder_roary, "gene_presence_absence.csv"), show_col_types = F) %>% clean_names()
gpa <- read_table(paste0(folder_roary, "gene_presence_absence.Rtab"), show_col_types = F)
gnew <- read_table(paste0(folder_roary, "number_of_new_genes.Rtab"), col_names = F, show_col_types = F)
gconserved <- read_table(paste0(folder_roary, "number_of_conserved_genes.Rtab"), col_names = F, show_col_types = F)
gunique <- read_table(paste0(folder_roary, "number_of_unique_genes.Rtab"), col_names = F, show_col_types = F)
gpan <- read_table(paste0(folder_roary, "number_of_genes_in_pan_genome.Rtab"), col_names = F, show_col_types = F)

# Remove rhizobium and keep only Ensifer
# gpa <- gpa %>%
#     select(-c(annotated_g1, annotated_g7, annotated_g12, annotated_g4, annotated_g18))


# pivot longer
gpa[gpa == 0] <- NA
gpa <- clean_names(gpa)
gpa_long <- gpa %>%
    pivot_longer(cols = -gene, values_drop_na = T, names_prefix = "annotated_", names_to = "genome_id") %>%
    mutate(genome_id = factor(genome_id, paste0("g", 1:19))) %>%
    clean_names()

# 0. some numbers
# Size of the Core Genome; gfs stands for gene frequency spectrum
gpa_glist <- gpa_long %>%
    group_by(gene) %>%
    count(name = "n_isolates") %>%
    mutate(n_isolates = factor(n_isolates, 1:19))

gpa_gfs <- gpa_glist %>%
    group_by(n_isolates) %>%
    count(name = "n_genes")
gpa_gfs$n_genes[gpa_gfs$n_isolates==19] # 1331 core genes when -i=75 and -cd=99

# Pan Genome Size
length(unique(gpa_glist$gene)) # 30408

# Size of the Accessory Genome
gpa_gfs$n_genes

#  number of singleton
gpa_gfs$n_genes[1] # 16913


# 1. Number of gene shared in these isolates
#' This is the gene frequency spectrum G(k)
#' that correlates the number of orthologous genes groups (OGGs) containing genes from exactly k genomes
p <- gpa_gfs %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_col(aes(x = n_isolates, y = n_genes), color = "black", fill = "white") +
    #scale_x_continuous(breaks = 1:8) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 18), minor_breaks = 1:18) +
    theme_bw() +
    theme(
    ) +
    guides() +
    labs(x = "gene shared by # of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-01-histogram.png"), plot = p, width = 5, height = 5)

# 2. heatmap for genes
gpa_glist_core <- gpa_glist %>% filter(n_isolates == 14)

p <- gpa_long %>%
    #filter(str_detect(gene, "b")) %>%
    mutate(value = factor(value)) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = gene, fill = value)) +
    scale_fill_manual(values = c(`1` = "maroon")) +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 10)
    ) +
    guides(fill = "none") +
    labs(x = "strain", y = "gene")
ggsave(paste0(folder_data, "temp/42-02-heatmap.png"), plot = p, width = 10, height = 20)


# 3. conserved gene vs. total genes in the pan genome
glist <- bind_rows(
    mutate(gconserved, key = "conserved genes", replicate = paste0("c", 1:10)),
    mutate(gpan, key = "total genes", replicate = paste0("t", 1:10))
) %>%
    pivot_longer(cols = -c(key, replicate), names_prefix = "X", names_to = "genome", values_to = "n_genes") %>%
    mutate(genome = as.numeric(genome), replicate = factor(replicate), key = factor(key, c("total genes", "conserved genes")))
p <- glist %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_line(aes(x = genome, y = n_genes, color = key, group = replicate), linewidth = .5) +
    scale_color_simpsons() +
    scale_x_continuous(minor_breaks = 1:20, limits = c(1,20)) +
    #scale_y_continuous(minor_breaks = 10, limits = c(1,20)) +
    theme_bw() +
    theme(
    ) +
    guides(color = guide_legend(title = "")) +
    labs(x = "# of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-03-conserved_vs_total_genes.png"), plot = p, width = 5, height = 3)

# numbers
range(glist$n_genes) # 1331 and 30408

# 4. new genes vs. unique genes
glist <- bind_rows(
    mutate(gunique, key = "unique genes", replicate = paste0("c", 1:10)),
    mutate(gnew, key = "new genes", replicate = paste0("t", 1:10))
) %>%
    pivot_longer(cols = -c(key, replicate), names_prefix = "X", names_to = "genome", values_to = "n_genes") %>%
    mutate(genome = as.numeric(genome), replicate = factor(replicate), key = factor(key, c("unique genes", "new genes")))
p <- glist %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_line(aes(x = genome, y = n_genes, color = key, group = replicate), linewidth = .5) +
    scale_color_simpsons() +
    scale_x_continuous(minor_breaks = 1:20, limits = c(1,20)) +
    #scale_y_continuous(minor_breaks = 10, limits = c(1,20)) +
    theme_bw() +
    theme(
    ) +
    guides(color = guide_legend(title = "")) +
    labs(x = "# of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-04-unique_vs_new_genes.png"), plot = p, width = 5, height = 3)

# Number
range(glist$n_genes) # 17 and 17166

# 5. all possible genes
glist <- bind_rows(
    mutate(gconserved, key = "conserved genes", replicate = paste0("c", 1:10)),
    mutate(gunique, key = "unique genes", replicate = paste0("u", 1:10)),
    mutate(gnew, key = "new genes", replicate = paste0("n", 1:10)),
    mutate(gpan, key = "total genes", replicate = paste0("t", 1:10)),
) %>%
    pivot_longer(cols = -c(key, replicate), names_prefix = "X", names_to = "genome", values_to = "n_genes") %>%
    mutate(genome = as.numeric(genome), replicate = factor(replicate), key = factor(key, c("total genes", "unique genes", "conserved genes", "new genes")))
p <- glist %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_line(aes(x = genome, y = n_genes, color = key, group = replicate), linewidth = .5) +
    scale_color_simpsons() +
    scale_x_continuous(minor_breaks = 1:20, limits = c(1,20)) +
    theme_bw() +
    theme(
    ) +
    guides(color = guide_legend(title = "")) +
    labs(x = "# of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-05-all_genes.png"), plot = p, width = 5, height = 3)

# 6. build a tree using the gene presence and absence data
library(vegan)
gpa_bin <- gpa
gpa_bin[is.na(gpa_bin)] <- 0
t_gpa <- gpa_bin %>%
    rename_with(~str_remove(., "^annotated_"), everything()) %>%
    select(-gene) %>%
    as.matrix() %>%
    t()
dim(t_gpa) # 19x30408

# Calculate Bray-Curtis dissimilarity
bin_dis <- vegdist(t_gpa, method = "bray")
clus <- hclust(bin_dis)
#dendrogram <- as.dendrogram(clus)
#plot(dendrogram)

library(ggtree)
library(ape)
tree <- as.phylo(clus)
p <- tree %>%
    ggtree() +
    geom_tiplab()
ggsave(paste0(folder_data, "temp/42-06-tree_all_genes.png"), plot = p, width = 6, height = 3)




