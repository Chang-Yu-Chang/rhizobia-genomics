#' This script cleans the tables of genetic distance and join them

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

genomes <- read_csv(paste0(folder_data, "temp/00-genomes.csv"))
folder_genomes <- paste0(folder_genomics, "genomes/")
genomes <- bind_rows(genomes, tibble(batch_name = rep("ncbi", 4), genome_name = c("em1021", "em1022", "usda1106", "wsm419"), genome_id = c("em1021", "em1022", "usda1106", "wsm419"), exp_id = c("em1021", "em1022", "usda1106", "wsm419")))

# 1. ANI 
dist_ani <- read_delim(paste0(folder_data, "genomics/popgen/fastani/ani.txt"), delim = "\t", col_names = F)
list_genomes <- read_delim(paste0(folder_data, "genomics/popgen/fastani/list_genomes.txt"), delim = "\t", col_names = F)

# 1.1 clean the tables
colnames(dist_ani) <- c("genome_id1", "genome_id2", "d_ani", "frag1", "frag2")
dist_ani <- dist_ani %>%
    # Compute so that d_ani represents the genomic distance instance of similarity
    mutate(d_ani = 1-d_ani/100) %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = ordered(genome_id1, genomes$genome_id), genome_id2 = ordered(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2) %>%
    filter(genome_id1 <= genome_id2) 
write_csv(dist_ani, paste0(folder_data, "temp/19-dist_ani.csv"))

# 2. k-mers  
# 2.1 genome kmers
dist_kmer <- read_delim(paste0(folder_data, "genomics/popgen/kmer/kmer.txt"))
list_sigs <- read_delim(paste0(folder_data, "genomics/popgen/kmer/list_sigs.txt"), delim = "\t", col_names = F)

dist_kmer <- dist_kmer %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_kmer")  %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = ordered(genome_id1, genomes$genome_id), genome_id2 = ordered(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2) %>%
    filter(genome_id1 <= genome_id2)

write_csv(dist_kmer, paste0(folder_data, "temp/19-dist_kmer.csv"))


# 2.2 contig kmers
dist_kmer_contigs <- read_delim(paste0(folder_data, "genomics/popgen/kmer_contigs/kmer.txt"))
list_sigs <- read_delim(paste0(folder_data, "genomics/popgen/kmer_contigs/list_sigs.txt"), delim = "\t", col_names = F)
list_sigs <- str_extract(list_sigs$X1, "g\\d+_contig_\\d+")
length(list_sigs) # 228 contigs
dist_kmer_contigs <- dist_kmer_contigs %>%
    mutate(contig_id1 = colnames(.)) %>%
    pivot_longer(-contig_id1, names_to = "contig_id2", values_to = "d_kmer")  %>%
    mutate(contig_id1 = str_remove(contig_id1, ".+/contigs/") %>% str_remove(".fasta")) %>%
    mutate(contig_id2 = str_remove(contig_id2, ".+/contigs/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = str_remove(contig_id1, "_contig_\\d"), genome_id2 = str_remove(contig_id2, "_contig_\\d")) %>%
    mutate(contig_id1 = ordered(contig_id1, list_sigs), contig_id2 = ordered(contig_id2, list_sigs)) %>%
    arrange(contig_id1, contig_id2) %>%
    filter(contig_id1 <= contig_id2) %>%
    select(starts_with("genome_id"), starts_with("contig"), d_kmer)

nrow(dist_kmer_contigs) # choose(228,2) + 228 = 26106

# 3. gene content
dist_jaccard <- read_csv(paste0(folder_data, "temp/13-dist_jaccard.csv"))
dist_fluidity <- read_csv(paste0(folder_data, "temp/13-dist_fluidity.csv"))
dist_jaccard_contigs <- read_csv(paste0(folder_data, "temp/13-dist_jaccard_contigs.csv"))

# 4. 825 core genes
dist_scccg <- read_csv(paste0(folder_data, "temp/15-dist_sccg.csv"))

# Join the distance tables
dist_genetics <- dist_ani %>%
    select(genome_id1, genome_id2, d_ani) %>%
    left_join(dist_kmer) %>%
    left_join(dist_jaccard) %>%    
    left_join(dist_fluidity) %>%
    left_join(dist_scccg)

write_csv(dist_genetics, paste0(folder_data, "temp/19-dist_genetics.csv"))


dist_genetics_contigs <- dist_kmer_contigs %>%
    left_join(dist_jaccard_contigs)

write_csv(dist_genetics_contigs, paste0(folder_data, "temp/19-dist_genetics_contigs.csv"))





if (FALSE) {


dist_jac <- read_csv(paste0(folder_data, "temp/17-dist_jac.csv"), show_col_types = F)
dist_flu <- read_csv(paste0(folder_data, "temp/17-dist_flu.csv"), show_col_types = F)

# 0.3 join 
dists <- dist_kmer %>%
    left_join(dist_jac) %>%
    left_join(dist_flu) %>%
    filter(genome_id1 < genome_id2)
write_csv(dists, paste0(folder_data, "temp/19-dists_genetic.csv"))


# 1. genetic distance 
p1 <- dists %>%
    ggplot() +
    geom_point(aes(x = d_kmer, y = d_jaccard), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    labs()

p2 <- dists %>%
    ggplot() +
    geom_point(aes(x = d_kmer, y = d_fluidity), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    labs()

p3 <- dists %>%
    ggplot() +
    geom_point(aes(x = d_jaccard, y = d_fluidity), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    labs()

p <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb")
ggsave(paste0(folder_data, "temp/19-01-g_dist.png"), p, width = 12, height = 4)

# 2. genetic distance colored by pair types 
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
    geom_point(aes(x = d_kmer, y = d_jaccard, color = pair_sp), shape = 21) +
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
    geom_point(aes(x = d_kmer, y = d_fluidity, color = pair_sp), shape = 21) +
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
    geom_point(aes(x = d_jaccard, y = d_fluidity, color = pair_sp), shape = 21) +
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

# 3. genetic distance colored by population 
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
    geom_point(aes(x = d_kmer, y = d_jaccard, color = pair_sp), shape = 21) +
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
    geom_point(aes(x = d_kmer, y = d_fluidity, color = pair_sp), shape = 21) +
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
    geom_point(aes(x = d_jaccard, y = d_fluidity, color = pair_sp), shape = 21) +
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
