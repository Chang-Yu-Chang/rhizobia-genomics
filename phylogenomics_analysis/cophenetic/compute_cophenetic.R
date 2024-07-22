#' This script plots the NJ trees at the contig level

renv::load()
library(tidyverse)
library(ggtree)
library(tidytree)
library(ape)
source(here::here("metadata.R"))

load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))

# Functions
get_group <- function (dm) {
    #' Get the group labels from dm
    tibble(genome_id = colnames(dm)) %>%
        left_join(isolates) %>%
        pull(site_group)
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


if (FALSE) {

# Example
set.seed(2)
tr <- rtree(10)
tr <- tr_gpa
plot(tr)
dm <- cophenetic(tr)
groups <- c(rep(c(1,2), 16))
obs <- calculate_distances(dm, groups) %>% calculate_dist_mean()
obs
pm <- permute_dist(dm, groups)
hist(pm)
}

# Master table
tb <- tibble(
    id = factor(1:4),
    feature = rep(c("core", "gene presence/absence"), each = 2),
    population = rep(c("VA", "PA"), 2),
    tree = list(
        drop.tip(tr, isolates$genome_id[isolates$population != "VA"]),
        drop.tip(tr, isolates$genome_id[isolates$population != "PA"]),
        drop.tip(tr_gpa, isolates$genome_id[isolates$population != "VA"]),
        drop.tip(tr_gpa, isolates$genome_id[isolates$population != "PA"])
    )
)

tb <- tb %>%
    rowwise() %>%
    mutate(
        dm = list(cophenetic(tree)),
        groups = list(get_group(dm)),
        distances_obs = list(calculate_distances(dm, groups) %>% calculate_dist_mean),
        distances_pm = list(permute_dist(dm, groups))
    )

tb_obs <- tb %>% unnest(distances_obs)

tb %>%
    #unnest(distances_obs) %>%
    unnest(distances_pm) %>%
    ggplot() +
    geom_histogram(aes(x = distances_pm), color = "black", fill = "white") +
    geom_vline(data = tb_obs, aes(xintercept = distances_obs), color = 2) +
    facet_grid(feature ~ population, scales = "free_x") +
    theme_bw() +
    theme() +
    guides() +
    labs()


compute_percentiles <- function(tb) {
    tb %>%
        group_by(id, feature, population) %>%
        arrange(distances_pm) %>%
        summarize(p05 = quantile(distances_pm, 0.05), p50 = quantile(distances_pm, 0.5), p95 = quantile(distances_pm, 0.95))
        #pivot_longer(-c(id, feature, population), names_to = "percentile")
}


p <- tb %>%
    #select(id, distances_pm) %>%
    unnest(distances_pm) %>%
    compute_percentiles() %>%
    ggplot() +
    geom_hline(yintercept = c(0, 1), color = "grey80", linetype = 2, linewidth = 1) +
    geom_segment(aes(x = id, xend = id, y = p05, yend = p95), linewidth = 1,  arrow = arrow(length = unit(3, "mm"), angle = 90, ends = "both")) +
    geom_point(aes(x = id, y = p50, color = "median"), shape = 21, stroke = 2, size = 2) +
    geom_point(data = tb_obs, aes(x = id, y = distances_obs, color = "obs"), shape = 21, stroke = 2, size = 2) +
    scale_x_discrete(breaks = 1:4, labels = tb$feature, position = "top") +
    scale_y_continuous(limits = c(0, 1.2)) +
    scale_color_manual(values = c("obs" = "maroon", "median" = "black"), name = NULL) +
    facet_grid(population~., switch = "y", scales = "free_y") +
    coord_flip() +
    theme_bw() +
    theme(
        legend.position = "top",
    ) +
    guides() +
    labs(x = "", y = "")

ggsave(paste0(folder_data, "phylogenomics_analysis/cophenetic/01-cophenetic.png"), p, width = 6, height = 4)






