#' This script make the distance matrix for NJ trees

renv::load()
library(tidyverse)
library(janitor)
library(ape)
library(tidytree)
source(here::here("analysis/00-metadata.R"))

#gd <- read_csv(paste0(folder_data, "genomics/pangenome/panaroo/gene_data.csv"))
dists <- read_csv(paste0(folder_data, 'temp/31-dists.csv'))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))
isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
list_traits <- c("dry_weight_mg", "nodule_number", "root_weight_mg",
        str_subset(colnames(isolates_traits), "r_"),
        str_subset(colnames(isolates_traits), "lag_"),
        str_subset(colnames(isolates_traits), "maxOD_"))
list_dists <- c("ani", "kmer", "jaccard", "fluidity", "hamming", "jc", "geo", "growth")

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

dists_species <- dists %>%
    mutate(species_pair = case_when(
        species1 == "medicae" & species2 == "medicae" ~ "medicae pair",
        species1 == "meliloti" & species2 == "meliloti" ~ "meliloti pair"
    ))

# Make NJ trees
make_tree <- function(di, d_trait) {
    #' This functions creates a phylo object from a long-format distance matrix
    if (d_trait %in% paste0("d_", list_traits)) {
        list_avails <- isolates_traits$genome_id[!is.na(isolates_traits[str_remove(d_trait, "d_")])]
        di <- bind_cols(select(di, genome_id1, genome_id2), tibble(dd = unlist(di[,d_trait]))) %>%
            filter(genome_id1 %in% list_avails, genome_id2 %in% list_avails)
    } else {
        di <- bind_cols(select(di, genome_id1, genome_id2), tibble(dd = unlist(di[,d_trait])))
    }
    colnames(di)[3] <- "dd"
    dist1 <- tibble(genome_id1 = di$genome_id1, genome_id2 = di$genome_id2, dd = di$dd)
    dist1_swapped <- tibble(genome_id1 = di$genome_id2, genome_id2 = di$genome_id1, dd = di$dd)
    dist1_swapped <- filter(dist1_swapped, genome_id1 != genome_id2)
    dist2 <- bind_rows(dist1, dist1_swapped) %>%
        drop_na(dd) %>%
        pivot_wider(names_from = genome_id2, values_from = dd) %>%
        select(-genome_id1) %>%
        as.matrix()
    tree <- dist2 %>%
        nj()
    return(tree)
}

# All 41 genomes
list_trees <- rep(list(NA), length(list_dists))
names(list_trees) <- list_dists
for (i in 1:length(list_dists)) list_trees[[i]] <- dists %>% make_tree(paste0("d_", list_dists[i]))

# meliloti
list_trees_meliloti <- rep(list(NA), length(list_dists))
names(list_trees_meliloti) <- list_dists
for (i in 1:length(list_dists)) {
    list_trees_meliloti[[i]] <- dists_species %>%
        filter(species_pair == "meliloti pair") %>%
        make_tree(paste0("d_", list_dists[i]))
}

# medicae
list_trees_medicae <- rep(list(NA), length(list_dists))
names(list_trees_medicae) <- list_dists
for (i in 1:length(list_dists)) {
    list_trees_medicae[[i]] <- dists_species %>%
        filter(species_pair == "medicae pair") %>%
        make_tree(paste0("d_", list_dists[i]))
}


# Make phylo object
write.tree(list_trees$jaccard, paste0(folder_data, "temp/32-jaccard.tre"))


# Make contig trees
dist_genetics_contigs <- read_csv(paste0(folder_data, 'temp/19-dist_genetics_contigs.csv'))
contigs <- read_csv(paste0(folder_data, 'temp/14-contigs.csv'))
dist_genetics_contigs <- dist_genetics_contigs %>%
    left_join(rename(rename_all(contigs, ~ paste0(.x, "1")))) %>%
    left_join(rename(rename_all(contigs, ~ paste0(.x, "2"))))


make_tree_c <- function(di, d_trait) {
    #' This functions creates a phylo object from a long-format distance matrix
    di <- bind_cols(select(di, contig_id1, contig_id2), tibble(dd = unlist(di[,d_trait])))
    colnames(di)[3] <- "dd"
    dist1 <- tibble(contig_id1 = di$contig_id1, contig_id2 = di$contig_id2, dd = di$dd)
    dist1_swapped <- tibble(contig_id1 = di$contig_id2, contig_id2 = di$contig_id1, dd = di$dd)
    dist1_swapped <- filter(dist1_swapped, contig_id1 != contig_id2)
    dist2 <- bind_rows(dist1, dist1_swapped) %>%
        drop_na(dd) %>%
        pivot_wider(names_from = contig_id2, values_from = dd) %>%
        select(-contig_id1) %>%
        as.matrix()
    tree <- dist2 %>%
        nj()
    return(tree)
}


list_trees_contigs <- list(
    ensifer_kmer = dist_genetics_contigs %>% make_tree_c("d_kmer"),
    meliloti_kmer =  dist_genetics_contigs %>% filter(species1 == "meliloti", species2 == "meliloti") %>% make_tree_c("d_kmer"),
    medicae_kmer =  dist_genetics_contigs %>% filter(species1 == "medicae", species2 == "medicae") %>% make_tree_c("d_kmer"),
    ensifer_jaccard = dist_genetics_contigs %>% make_tree_c("d_jaccard"),
    meliloti_jaccard =  dist_genetics_contigs %>% filter(species1 == "meliloti", species2 == "meliloti") %>% make_tree_c("d_jaccard"),
    medicae_jaccard =  dist_genetics_contigs %>% filter(species1 == "medicae", species2 == "medicae") %>% make_tree_c("d_jaccard")
)


save(list_trees, list_trees_meliloti, list_trees_medicae, list_trees_contigs,
    file = paste0(folder_data, "temp/32-trees.RData"))

if (FALSE) {
    library(TreeDist) # For computing RF distance
    # Compute the incongruence among trees using RF distance
    TreeDistance(list_trees, list_trees)
    TreeDistance(list_trees_meliloti, list_trees_meliloti)
    TreeDistance(list_trees_medicae, list_trees_medicae)

}
