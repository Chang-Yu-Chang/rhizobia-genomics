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

#
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
        scale_x_continuous(breaks = c(0,0.5,1)) +
        scale_y_continuous(breaks = c(0,0.5,1)) +
        facet_grid(~d_type) +
        theme_bw() +
        labs(x = "d_genomic", y = d_trait)
}

p1 <- plot_dots(dists_long, "d_growth") +
theme(axis.title.x = element_blank(), strip.text.x = element_text(size = 15), strip.background = element_rect(color = NA, fill = "white")) +
    labs(x = NULL)
p2 <- plot_dots(dists_long, "d_symbiosis") + 
    theme(strip.text = element_blank(),plot.margin = unit(c(0,0,0,0), "mm"))
    
p <- plot_grid(p1, p2, nrow = 2, align = "vh", axis = "lrtb", scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-01-genetics_vs_phenotypes.png"), p, width = 12, height = 6)

# 2. genetic distance vs one trait
p1 <- plot_dots(dists_long, "d_r_30c") +
    theme(axis.title.x = element_blank(), strip.text.x = element_text(size = 15), strip.background = element_rect(color = NA, fill = "white")) +
    labs(x = NULL)
p2 <- plot_dots(dists_long, "d_lag_30c") +
    theme(strip.text = element_blank(),plot.margin = unit(c(0,0,0,0), "mm")) +
    labs(x = NULL)
p3 <- plot_dots(dists_long, "d_nodule_number") +
    theme(strip.text = element_blank(),plot.margin = unit(c(0,0,0,0), "mm")) +
    labs(x = NULL)
p4 <- plot_dots(dists_long, "d_dry_weight_mg") +
    theme(strip.text = element_blank(),plot.margin = unit(c(0,0,0,0), "mm"))
p <- plot_grid(p1, p2, p3, p4,
    ncol = 1, align = "v", axis = "lr", scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-02-genetics_vs_trait.png"), p, width = 12, height = 12)

# 3. Correlation between traits
list_genomics <- c("ani", "kmer", "jaccard", "fluidity")
list_traits <- c("dry_weight_mg", "nodule_number", "root_weight_mg",
        str_subset(colnames(isolates_traits), "r_"),
        str_subset(colnames(isolates_traits), "lag_"),
        str_subset(colnames(isolates_traits), "maxOD_"))

tb_gt <- crossing(d_g = list_genomics, d_t = list_traits) %>%
    rowwise %>%
    mutate(dat = list(tibble(dists[,paste0("d_", c(d_g, d_t))]))) %>%
    ungroup()


tb_gt %>%
    map(~cor(.x$dat))




cor.test(dists$d_kmer, dists$d_growth)
cor.test(dists$d_kmer, dists$d_symbiosis)
cor.test(dists$d_kmer, dists$d_dry_weight_mg)
cor.test(dists$d_kmer, dists$d_nodule_number)
