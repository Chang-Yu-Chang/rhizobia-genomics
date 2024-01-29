#' This script cleans the tables of genetic distance and join them

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

genomes <- read_csv(paste0(folder_data, "temp/00-genomes.csv"))
folder_genomes <- paste0(folder_genomics, "genomes/")

# 1. ANI 
dist_ani <- read_delim(paste0(folder_data, "genomics/popgen/fastani/ani.txt"), delim = "\t", col_names = F)
list_genomes <- read_delim(paste0(folder_data, "genomics/popgen/fastani/list_genomes.txt"), delim = "\t", col_names = F)

# 1.1 clean the tabkes
colnames(dist_ani) <- c("genome_id1", "genome_id2", "distance_ani", "frag1", "frag2")
dist_ani <- dist_ani %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = factor(genome_id1, genomes$genome_id), genome_id2 = factor(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2)

write_csv(dist_ani, paste0(folder_data, "temp/13-dist_ani.csv"))

# 2. k-mers 
dist_kmer <- read_delim(paste0(folder_data, "genomics/popgen/kmer/kmer.txt"))
list_sigs <- read_delim(paste0(folder_data, "genomics/popgen/kmer/list_sigs.txt"), delim = "\t", col_names = F)

# 2.1 clean the tables
dist_kmer <- dist_kmer %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "distance_kmer")  %>%
        mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = factor(genome_id1, genomes$genome_id), genome_id2 = factor(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2)

write_csv(dist_kmer, paste0(folder_data, "temp/13-dist_kmer.csv"))

# 3. genomic fluidity

if (FALSE) {


dist_jac <- read_csv(paste0(folder_data, "temp/17-dist_jac.csv"), show_col_types = F)
dist_flu <- read_csv(paste0(folder_data, "temp/17-dist_flu.csv"), show_col_types = F)

# 0.3 join ----
dists <- dist_kmer %>%
    left_join(dist_jac) %>%
    left_join(dist_flu) %>%
    filter(genome_id1 < genome_id2)
write_csv(dists, paste0(folder_data, "temp/19-dists_genetic.csv"))


# 1. genetic distance ----
p1 <- dists %>%
    ggplot() +
    geom_point(aes(x = distance_kmer, y = distance_jaccard), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    labs()

p2 <- dists %>%
    ggplot() +
    geom_point(aes(x = distance_kmer, y = distance_fluidity), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    labs()

p3 <- dists %>%
    ggplot() +
    geom_point(aes(x = distance_jaccard, y = distance_fluidity), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    labs()

p <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb")
ggsave(paste0(folder_data, "temp/19-01-g_dist.png"), p, width = 12, height = 4)

# 2. genetic distance colored by pair types ----
dists_meta <- dists %>%
    left_join(rename_with(isolates_mash, function (x) paste0(x, 1))) %>%
    left_join(rename_with(isolates_mash, function (x) paste0(x, 2)))

dists_meta <- dists_meta %>%
    mutate(pair_sp = case_when(
        species_name1 == "Ensifer meliloti" & species_name2 == "Ensifer meliloti" ~ "conspecific",
        species_name1 == "Ensifer medicae" & species_name2 == "Ensifer medicae" ~ "conspecific",
        species_name1 == "Ensifer sp." & species_name2 == "Ensifer sp." ~ "conspecific",
        species_name1 == "Ensifer adhaerens" & species_name2 == "Ensifer adhaerens" ~ "conspecific",

        species_name1 == "Ensifer meliloti" & species_name2 == "Ensifer medicae" ~ "meliloti & medicae",
        species_name1 == "Ensifer medicae" & species_name2 == "Ensifer meliloti" ~ "meliloti & medicae",
        TRUE ~ "others"
        # species_name1 == "Ensifer meliloti" & species_name2 == "Ensifer sp." ~ "meliloti & sp",
        # species_name1 == "Ensifer sp." & species_name2 == "Ensifer meliloti" ~ "meliloti & sp",
        # species_name1 == "Ensifer meliloti" & species_name2 == "Ensifer adhaerens" ~ "meliloti & sp",
        # species_name1 == "Ensifer adhaerens" & species_name2 == "Ensifer meliloti" ~ "meliloti & sp",
        #
        #
        # species_name1 == "Ensifer medicae" & species_name2 == "Ensifer sp." ~ "medicae & sp",
        # species_name1 == "Ensifer sp." & species_name2 == "Ensifer medicae" ~ "medicae & sp",
        # species_name1 == "Ensifer medicae" & species_name2 == "Ensifer adhaerens" ~ "medicae & sp",
        # species_name1 == "Ensifer adhaerens" & species_name2 == "Ensifer medicae" ~ "medicae & sp"
    ))

p1 <- dists_meta %>%
    ggplot() +
    geom_point(aes(x = distance_kmer, y = distance_jaccard, color = pair_sp), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_color_npg() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
    ) +
    guides(color = "none") +
    labs()
p2 <- dists_meta %>%
    ggplot() +
    geom_point(aes(x = distance_kmer, y = distance_fluidity, color = pair_sp), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_color_npg() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs()
p3 <- dists_meta %>%
    ggplot() +
    geom_point(aes(x = distance_jaccard, y = distance_fluidity, color = pair_sp), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_color_npg() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides(color = "none") +
    labs()

p <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb")
ggsave(paste0(folder_data, "temp/19-02-g_dist_sp.png"), p, width = 12, height = 4)

# 3. genetic distance colored by population ----
dists_meta <- dists %>%
    left_join(rename_with(isolates_mash, function (x) paste0(x, 1))) %>%
    left_join(rename_with(isolates_mash, function (x) paste0(x, 2)))

dists_meta <- dists_meta %>%
    mutate(pair_sp = case_when(
        rhizobia_population1 == "MLBS" & rhizobia_population2 == "MLBS" ~ "within MLBS",
        rhizobia_population1 == "Phila" & rhizobia_population2 == "Phila" ~ "within Phila",
        rhizobia_population1 == "MLBS" & rhizobia_population2 == "Phila" ~ "between MLBS and Phila",
        rhizobia_population1 == "Phila" & rhizobia_population2 == "MLBS" ~ "between MLBS and Phila",
        TRUE ~ "others"
    ))

p1 <- dists_meta %>%
    ggplot() +
    geom_point(aes(x = distance_kmer, y = distance_jaccard, color = pair_sp), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_color_npg() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides(color = "none") +
    labs()
p2 <- dists_meta %>%
    ggplot() +
    geom_point(aes(x = distance_kmer, y = distance_fluidity, color = pair_sp), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_color_npg() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs()
p3 <- dists_meta %>%
    ggplot() +
    geom_point(aes(x = distance_jaccard, y = distance_fluidity, color = pair_sp), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_color_npg() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides(color = "none") +
    labs()

p <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb")
ggsave(paste0(folder_data, "temp/19-03-g_dist_pop.png"), p, width = 12, height = 4)

} 
