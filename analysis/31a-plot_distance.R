#' This script plots the distances

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
library(ggsci)
source(here::here("analysis/00-metadata.R"))

isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))
dists <- read_csv(paste0(folder_data, "temp/31-dists.csv"))
dists_long <- read_csv(paste0(folder_data, "temp/31-dists_long.csv"))

nrow(dists) # choose(42,2)+42  = 903

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

p1 <- dists_long %>%
    filter(d_type != "fluidity") %>%
    plot_dots("d_growth") +
    theme(axis.title.x = element_blank(), strip.text.x = element_text(size = 13), strip.background = element_rect(color = NA, fill = "white")) +
    labs(x = NULL)
p2 <- dists_long %>%
    filter(d_type != "fluidity") %>%
    plot_dots("d_symbiosis") + 
    theme(axis.title.x = element_blank(), strip.text.x = element_text(size = 13), strip.background = element_rect(color = NA, fill = "white"))
p <- plot_grid(p1, p2, nrow = 2, align = "vh", axis = "lrtb", scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-01-genetics_vs_phenotypes.png"), p, width = 15, height = 6)

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


# 3. plot the mean value of two composite traits
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

# 4. plot the pairwise comparison between genetic distance s
plot_two_dists <- function (dists, dist1, dist2) {
    dist_i <- dists[,c(dist1, dist2)] %>% setNames(c("d1", "d2"))
    dist_i %>%
        ggplot() +
        geom_point(aes(x = d1, y = d2), shape = 1, size = 2, alpha = 0.5) +
        geom_smooth(aes(x = d1, y = d2), method = "lm") +
        scale_x_continuous(breaks = c(0,0.5,1)) +
        scale_y_continuous(breaks = c(0,0.5,1)) +
        theme_bw() +
        theme() +
        labs(x = dist1, y = dist2)
}

p1 <- dists %>% plot_two_dists("d_kmer", "d_ani")
p2 <- dists %>% plot_two_dists("d_jaccard", "d_ani")
p3 <- dists %>% plot_two_dists("d_hamming", "d_ani")
p4 <- dists %>% plot_two_dists("d_jc", "d_ani")
p5 <- dists %>% plot_two_dists("d_jaccard", "d_kmer")
p6 <- dists %>% plot_two_dists("d_hamming", "d_kmer")
p7 <- dists %>% plot_two_dists("d_jc", "d_kmer")
p8 <- dists %>% plot_two_dists("d_hamming", "d_jaccard")
p9 <- dists %>% plot_two_dists("d_jc", "d_jaccard")
p10 <- dists %>% plot_two_dists("d_jc", "d_hamming")


p <- plot_grid(p1, p2, p3, p4, NULL, p5, p6, p7, NULL, NULL, p8, p9, NULL, NULL, NULL, p10,
    nrow = 4, align = "hv", axis = "tblr", scale = 0.95, labels = c("A", "B", "C", "D", " ", "E", "F", "G", "", "", "H", "I", "", "", "", "J")
) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-04-genetic_dist_pairs.png"), p, width = 12, height = 12)

# 5. Plot the genetic distances vs phenotypic distances, colored by population
isolates_for_pair <- left_join(
    isolates_contigs %>% select(genome_id, species),
    isolates_traits %>% select(genome_id, population)
) %>%
    filter(!genome_id %in% c("g20", "g28"))

dists <- dists %>%
    left_join(rename(isolates_for_pair, genome_id1 = genome_id, species1 = species, population1 = population)) %>%
    left_join(rename(isolates_for_pair, genome_id2 = genome_id, species2 = species, population2 = population)) %>%
    mutate(compare_species = ifelse(species1 == species2, "within species", "between species")) %>%
    mutate(compare_population = ifelse(population1 == population2, "within population", "between populations"))

plot_two_dists_color <- function (dists, dist1, dist2, by_color) {
    dist_i <- dists[,c(dist1, dist2, by_color)] %>% setNames(c("d1", "d2", "bc"))
    dist_i %>%
        drop_na() %>%
        ggplot() +
        geom_point(aes(x = d1, y = d2, color = bc), shape = 1, size = 2, alpha = 0.5) +
        #geom_smooth(aes(x = d1, y = d2, color = bc), method = "lm") +
        scale_x_continuous(breaks = c(0,0.5,1)) +
        scale_y_continuous(breaks = c(0,0.5,1)) +
        scale_color_npg() +
        theme_bw() +
        theme(
            legend.position = "top"
        ) +
        guides(color = guide_legend(title = NULL)) +
        labs(x = dist1, y = dist2)
}

p1 <- dists %>% plot_two_dists_color("d_ani", "d_growth", "compare_population")
p2 <- dists %>% plot_two_dists_color("d_kmer", "d_growth", "compare_population")
p3 <- dists %>% plot_two_dists_color("d_jaccard", "d_growth", "compare_population")
p4 <- dists %>% plot_two_dists_color("d_hamming", "d_growth", "compare_population")
p5 <- dists %>% plot_two_dists_color("d_jc", "d_growth", "compare_population")

p6 <- dists %>% plot_two_dists_color("d_ani", "d_symbiosis", "compare_population")
p7 <- dists %>% plot_two_dists_color("d_kmer", "d_symbiosis", "compare_population")
p8 <- dists %>% plot_two_dists_color("d_jaccard", "d_symbiosis", "compare_population")
p9 <- dists %>% plot_two_dists_color("d_hamming", "d_symbiosis", "compare_population")
p10 <- dists %>% plot_two_dists_color("d_jc", "d_symbiosis", "compare_population")

p <- plot_grid(p1, p2,p3,p4,p5,p6,p7,p8,p9,p10,
    nrow = 2, align = "vh", axis = "lrtb", scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-05-genetics_vs_phenotypes.png"), p, width = 15, height = 6)

# 6. Plot the genetic distances vs phenotypic distances, colored by species
p1 <- dists %>% plot_two_dists_color("d_ani", "d_growth", "compare_species")
p2 <- dists %>% plot_two_dists_color("d_kmer", "d_growth", "compare_species")
p3 <- dists %>% plot_two_dists_color("d_jaccard", "d_growth", "compare_species")
p4 <- dists %>% plot_two_dists_color("d_hamming", "d_growth", "compare_species")
p5 <- dists %>% plot_two_dists_color("d_jc", "d_growth", "compare_species")

p6 <- dists %>% plot_two_dists_color("d_ani", "d_symbiosis", "compare_species")
p7 <- dists %>% plot_two_dists_color("d_kmer", "d_symbiosis", "compare_species")
p8 <- dists %>% plot_two_dists_color("d_jaccard", "d_symbiosis", "compare_species")
p9 <- dists %>% plot_two_dists_color("d_hamming", "d_symbiosis", "compare_species")
p10 <- dists %>% plot_two_dists_color("d_jc", "d_symbiosis", "compare_species")

p <- plot_grid(p1, p2,p3,p4,p5,p6,p7,p8,p9,p10,
    nrow = 2, align = "vh", axis = "lrtb", scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/31a-06-genetics_vs_phenotypes.png"), p, width = 15, height = 6)
