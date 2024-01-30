#' This script plots the distances

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
dists <- read_csv(paste0(folder_data, "temp/31-dists.csv"))
dists_long <- read_csv(paste0(folder_data, "temp/31-dists_long.csv"))

nrow(dists) # 903
choose(42,2)+42 

dists_long <- dists_long %>% 
    mutate(d_type = factor(d_type, unique(dists_long$d_type)))

# 1. genetic distances vs phenotypic distances
plot_dots <- function (dists_long, d_trait) {
    #' Plot the scatterplot of genomic distance agains a trait distance
    #' d_trait is a string
    dist_trait <- select(dists, genome_id1, genome_id2)
    dist_trait$dd <- unlist(dists[d_trait])
    dists_long %>%
        filter(d_group == "genetic") %>% 
        left_join(dist_trait) %>%
        drop_na() %>%
        ggplot() +
        geom_point(aes(x = value, y = dd), shape = 1, size = 2, alpha = 0.5) +
        geom_smooth(aes(x = value, y = dd), method = "lm") +
        scale_x_continuous(breaks = c(0,0.5,1)) +
        scale_y_continuous(breaks = c(0,0.5,1)) +
        facet_grid(~d_type) +
        theme_bw() +
        theme(axis.title.x = element_blank(), strip.text.x = element_text(size = 15), strip.background = element_rect(color = NA, fill = "white")) +
        labs(x = "d_genomic", y = d_trait)
}

p1 <- plot_dots(dists_long, "d_growth") +
    theme(axis.title.x = element_blank(), strip.text.x = element_text(size = 13), strip.background = element_rect(color = NA, fill = "white")) +
    labs(x = NULL)
p2 <- plot_dots(dists_long, "d_symbiosis") + 
    theme(plot.margin = unit(c(0,0,0,0), "mm")) 
p <- plot_grid(p1, p2, nrow = 2, align = "vh", axis = "lrtb", scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-01-genetics_vs_phenotypes.png"), p, width = 12, height = 6)

# 2. genetic distance vs one trait
p1 <- plot_dots(dists_long, "d_r_30c") +
    labs(x = NULL)
p2 <- plot_dots(dists_long, "d_lag_30c") +
    theme(plot.margin = unit(c(0,0,0,0), "mm")) +
    labs(x = NULL)
p3 <- plot_dots(dists_long, "d_r_25c") +
    theme(plot.margin = unit(c(0,0,0,0), "mm")) +
    labs(x = NULL)
p4 <- plot_dots(dists_long, "d_nodule_number") +
    theme(plot.margin = unit(c(0,0,0,0), "mm")) +
    labs(x = NULL)
p5 <- plot_dots(dists_long, "d_dry_weight_mg") +
    theme(plot.margin = unit(c(0,0,0,0), "mm")) +
    labs(x = NULL)
p6 <- plot_dots(dists_long, "d_root_weight_mg") +
    theme(plot.margin = unit(c(0,0,0,0), "mm"))

p <- plot_grid(p1, p2, p3, p4, p5, p6, 
    ncol = 1, align = "v", axis = "lr", scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-02-genetics_vs_trait.png"), p, width = 12, height = 18)


# 3. plot the two composite traits
dists_i <- dists %>%
    filter(genome_id1 != genome_id2) %>%
    select(d_growth, d_symbiosis) %>%
    drop_na()
p <- dists_i %>%
    ggplot() +
    geom_point(aes(x = d_growth, y = d_symbiosis), shape = 1, size = 2, alpha = 0.8) +
    geom_smooth(aes(x = d_growth, y = d_symbiosis), method = "lm") +
    scale_x_continuous(breaks = c(0,0.5,1)) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    theme_bw() +
    theme() +
    labs()
ggsave(paste0(folder_data, "temp/31a-03-growth_vs_symbiosis.png"), p, width = 4, height = 4)
cor.test(dists_i$d_growth, dists_i$d_symbiosis, method = "spearman")
