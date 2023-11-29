#' This script permute the genomes included in a pan genome analysis

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 0. read data
egc <- read_delim(paste0(folder_data, "/temp/anvio/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)

# Extract only the unique gene id
egc_g <- egc %>%
    select(unique_id, gene_cluster_id, genome_name) %>%
    distinct(gene_cluster_id, genome_name, .keep_all = T) %>%
    left_join(genomes) %>%
    select(-genome_name)

n_b = 100 # number of permutation
n_g = 17 # total number of genomes in the pangenome

list_pangenome <- rep(list(NA), n_g-1)

for (j in 1:(n_g-1)) { # number of genomes sampled
    list_boot <- rep(list(NA), n_b)
    for (i in 1:n_b) {
        set.seed(i)
        g_sampled <- sample(genomes$genome_id, size = j, replace = F)

        gfs <- egc_g %>% # gene frequency spectrum
            filter(genome_id %in% g_sampled) %>%
            group_by(gene_cluster_id) %>%
            summarize(n_genomes = n()) %>%
            group_by(n_genomes) %>%
            summarize(n_genes = n())

        g1 <- gfs$n_genes[gfs$n_genomes == 1] # singleton
        g2 <- gfs$n_genes[gfs$n_genomes == 2] # duplicate
        gc <- gfs$n_genes[gfs$n_genomes == j] # core gene clusters
        gp <- sum(gfs$n_genes) # total number of gene clusters in the pangenome

        list_boot[[i]] = tibble(bootstrap = i, singleton = g1, duplicate = g2, core = gc, pangene = gp)
    }

    list_pangenome[[j]]  <- bind_rows(list_boot) %>%
        mutate(n_genomes_sampled = j) %>%
        select(n_genomes_sampled, everything())

    print(j)
}

pangenome_boots <- bind_rows(list_pangenome)

write_csv(pangenome_boots, file = paste0(folder_data, "temp/44a-pangenome_boots.csv"))












