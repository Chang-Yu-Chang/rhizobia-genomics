#' This script plots the NJ trees at the contig level

renv::load()
library(tidyverse)
library(ggtree)
library(tidytree)
library(ape)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv")) %>%
    left_join(isolates) %>%
    select(genome_id, contig_id, species, replicon, replicon_type, population, site, site_group)
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))

# Functions
get_group <- function (dm, by_what = "species") {
    #' Get the group labels from dm
    tibble(genome_id = colnames(dm)) %>%
        mutate(genome_id = str_remove(genome_id, "_contig_\\d+")) %>%
        left_join(isolates_contigs) %>%
        left_join(isolates) %>%
        `[`(by_what) %>%
        unlist
}

calculate_distances <- function(dist_matrix, groups) {
    #' Compute the between cluster distance vs the within cluster distance
    n <- nrow(dist_matrix)
    within_distances <- c()
    between_distances <- c()

    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            if (groups[i] == groups[j]) {
                within_distances <- c(within_distances, dist_matrix[i, j])
            } else {
                between_distances <- c(between_distances, dist_matrix[i, j])
            }
        }
    }

    list(within = within_distances, between = between_distances)
}

calculate_dist_mean <- function(distances) {
    mean(distances$within) / mean(distances$between)
}

permute_dist <- function (dm, groups, np = 1000) {
    #' Permute the distance matrix and compute the within vs between cluster distance
    set.seed(123)  # For reproducibility
    distances <- calculate_distances(dm, groups)
    observed_diff <- mean(distances$within) - mean(distances$between)
    perm_diffs <- numeric(np)

    for (i in 1:np) {
        permuted_groups <- sample(groups)
        perm_distances <- calculate_distances(dm, permuted_groups)
        # Distance = within / between groups
        perm_diffs[i] <- mean(perm_distances$within) / mean(perm_distances$between)
    }
    return(perm_diffs)
}

compute_percentiles <- function(tb) {
    tb %>%
        group_by(id, feature, population) %>%
        arrange(distances_pm) %>%
        summarize(p05 = quantile(distances_pm, 0.05), p50 = quantile(distances_pm, 0.5), p95 = quantile(distances_pm, 0.95))
}

# Master table
tb <- tibble(
    id = factor(1:12),
    feature = rep(c("core", rep("gene presence/absence", 4), "structural variants"), each = 2),
    replicon_type = rep(c("core", "genome", "chromosome", "psyma", "psymb", "structural variants"), each = 2),
    population = rep(c("VA", "PA"), 6),
    tree = list(
        drop.tip(tr_seq_core, isolates$genome_id[isolates$population != "VA"]),
        drop.tip(tr_seq_core, isolates$genome_id[isolates$population != "PA"]),
        drop.tip(tr_gpa_genomes, isolates$genome_id[isolates$population != "VA"]),
        drop.tip(tr_gpa_genomes, isolates$genome_id[isolates$population != "PA"]),
        drop.tip(tr_gpa_chrom, contigs$contig_id[contigs$population != "VA"]),
        drop.tip(tr_gpa_chrom, contigs$contig_id[contigs$population != "PA"]),
        drop.tip(tr_gpa_psyma, contigs$contig_id[contigs$population != "VA"]),
        drop.tip(tr_gpa_psyma, contigs$contig_id[contigs$population != "PA"]),
        drop.tip(tr_gpa_psymb, contigs$contig_id[contigs$population != "VA"]),
        drop.tip(tr_gpa_psymb, contigs$contig_id[contigs$population != "PA"]),
        drop.tip(tr_spa_genomes, isolates$genome_id[isolates$population != "VA"]),
        drop.tip(tr_spa_genomes, isolates$genome_id[isolates$population != "PA"])
    )
)


# 1. grouping by species
tb1 <- tb %>%
    rowwise() %>%
    mutate(
        dm = list(cophenetic(tree)),
        groups = list(get_group(dm, by_what = "species")),
        distances_obs = list(calculate_distances(dm, groups) %>% calculate_dist_mean),
        distances_pm = list(permute_dist(dm, groups))
    )

tb_obs <- tb1 %>% unnest(distances_obs)

p <- tb1 %>%
    mutate(id = factor(id, 12:1)) %>%
    unnest(distances_pm) %>%
    compute_percentiles() %>%
    ggplot() +
    geom_hline(yintercept = c(0, 1), color = "grey80", linetype = 2, linewidth = 1) +
    geom_segment(aes(x = id, xend = id, y = p05, yend = p95), linewidth = 1,  arrow = arrow(length = unit(3, "mm"), angle = 90, ends = "both")) +
    geom_point(aes(x = id, y = p50, color = "median"), shape = 3, stroke = 1, size = 2) +
    geom_point(data = tb_obs, aes(x = id, y = distances_obs, color = "obs"), shape = 21, stroke = 2, size = 2) +
    scale_x_discrete(breaks = 1:12, labels = tb$replicon_type, position = "top") +
    scale_color_manual(values = c("obs" = "maroon", "median" = "black"), name = NULL) +
    facet_grid(population~., switch = "y", scales = "free_y") +
    coord_flip() +
    theme_bw() +
    theme(
        legend.position = "top",
    ) +
    guides() +
    labs(x = "", y = "meliloti vs medicae")

ggsave(paste0(folder_data, "phylogenomics_analysis/cophenetic/01-cophenetic_species.png"), p, width = 6, height = 4)

# 2. grouping by population
tb2 <- tb %>%
    rowwise() %>%
    mutate(
        dm = list(cophenetic(tree)),
        groups = list(get_group(dm, by_what = "species")),
        distances_obs = list(calculate_distances(dm, groups) %>% calculate_dist_mean),
        distances_pm = list(permute_dist(dm, groups))
    )

tb_obs <- tb2 %>% unnest(distances_obs)

p <- tb2 %>%
    mutate(id = factor(id, 12:1)) %>%
    unnest(distances_pm) %>%
    compute_percentiles() %>%
    ggplot() +
    geom_hline(yintercept = c(0, 1), color = "grey80", linetype = 2, linewidth = 1) +
    geom_segment(aes(x = id, xend = id, y = p05, yend = p95), linewidth = 1,  arrow = arrow(length = unit(3, "mm"), angle = 90, ends = "both")) +
    geom_point(aes(x = id, y = p50, color = "median"), shape = 3, stroke = 1, size = 2) +
    geom_point(data = tb_obs, aes(x = id, y = distances_obs, color = "obs"), shape = 21, stroke = 2, size = 2) +
    scale_x_discrete(breaks = 1:12, labels = tb$replicon_type, position = "top") +
    scale_color_manual(values = c("obs" = "maroon", "median" = "black"), name = NULL) +
    facet_grid(population~., switch = "y", scales = "free_y") +
    coord_flip() +
    theme_bw() +
    theme(
        legend.position = "top",
    ) +
    guides() +
    labs(x = "", y = "high vs low or urban vs suburban")

ggsave(paste0(folder_data, "phylogenomics_analysis/cophenetic/02-cophenetic_population.png"), p, width = 6, height = 4)






