#' This script joins the genetic distance

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    library(ggsci)
    #library(vcfR) # for handling VCF
    #library(poppr) # for pop gen analysis
    #library(ggtree)
    source(here::here("analysis/00-metadata.R"))
})

genome_kmer <- read_delim(paste0(folder_data, "genomics_old/popgen/genome_kmer/genome_kmer.txt"), show_col_types = F)
list_sig <- read_delim(paste0(folder_data, "genomics_old/popgen/genome_kmer/list_sig.txt"), delim = "\t", col_names = "file_name", show_col_types = F)
isolates_mash <- read_csv(paste0(folder_data, "temp/14-isolates_mash.csv"), show_col_types = F)

# 0. clean names ----
# 0.1 kmers ----
list_sig <- list_sig %>%
    mutate(genome_name = str_remove(file_name, ".*/genomes/") %>% str_remove("/04-taxonomy.*")) %>%
    left_join(isolates_mash)

genome_kmer <- 1-genome_kmer
genome_kmer <- genome_kmer %>%
    mutate(row_name = colnames(.)) %>%
    pivot_longer(-row_name, names_to = "col_name", values_to = "distance") %>%
    # Clean the matrix name
    mutate(row_name = str_remove(row_name, ".*/genomes/") %>% str_remove("/02.*")) %>%
    mutate(col_name = str_remove(col_name, ".*/genomes/") %>% str_remove("/02.*")) %>%
    pivot_wider(names_from = col_name, values_from = distance)

dist_kmer <- genome_kmer %>%
    pivot_longer(-row_name, names_to = "col_name", values_to = "distance_kmer") %>%
    left_join(rename_with(list_sig, function (x) paste0(x, 1)) %>% rename(col_name = genome_name1)) %>%
    left_join(rename_with(list_sig, function (x) paste0(x, 2)) %>% rename(row_name = genome_name2)) %>%
    select(genome_id1, genome_id2, row_name, col_name, distance_kmer)

# 0.2 jaccard distance ----
dist_jac <- read_csv(paste0(folder_data, "temp/17-dist_jac.csv"), show_col_types = F)

# 0.2 genomic fluidity ----
dist_flu <- read_csv(paste0(folder_data, "temp/17-dist_flu.csv"), show_col_types = F)

# 0.3 join ----
dists <- dist_kmer %>%
    left_join(dist_jac) %>%
    left_join(dist_flu) %>%
    filter(genome_id1 < genome_id2)

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
        rhizobia_population1 == "Phila" & rhizobia_population2 == "MLBS" ~ "conspecific",
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











