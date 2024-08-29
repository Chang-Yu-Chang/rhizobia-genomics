#' This script aggregates the pairwise genomic distance and perform mantel test

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(ape) # for eading trees
library(vegan) # for jaccard distance and mantel test
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))

# 1. single copy core genes ----
tr_seq_core <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/mltree/seq_core/seq_core.treefile"))
dist_sccg <- cophenetic.phylo(tr_seq_core) %>%
    as_tibble %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_sccg") %>%
    mutate(genome_id1 = factor(genome_id1, genomes$genome_id), genome_id2 = factor(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2)

# 2. gene content ----
## jaccard distance in gene presence and absence
gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
pat <- as.matrix(gpat[,-1])
dim(pat) # 36 x 26504
rownames(pat) <- gpat$genome_id
jdm <- vegdist(pat, method = "jaccard")

dist_jaccard <- jdm %>%
    as.matrix() %>% as_tibble() %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_jaccard")

# 3. ANI ----
dist_ani <- read_delim(paste0(folder_data, "genomics/ani/ani_genomes.txt"), delim = "\t", col_names = F)
list_genomes <- read_delim(paste0(folder_data, "genomics/ani/list_genomes.txt"), delim = "\t", col_names = F)

colnames(dist_ani) <- c("genome_id1", "genome_id2", "d_ani", "frag1", "frag2")
dist_ani <- dist_ani %>%
    # Compute so that d_ani represents the genomic distance instance of similarity
    mutate(d_ani = 1-d_ani/100) %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = ordered(genome_id1, genomes$genome_id), genome_id2 = ordered(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2) %>%
    select(genome_id1, genome_id2, d_ani)

# 4. k-mers ----
dist_kmer <- read_delim(paste0(folder_data, "genomics/kmer/kmer.txt"))
list_sigs <- read_delim(paste0(folder_data, "genomics/kmer/list_sigs.txt"), delim = "\t", col_names = F)

dist_kmer <- dist_kmer %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_kmer")  %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = ordered(genome_id1, genomes$genome_id), genome_id2 = ordered(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2)



# Join the distance tables ----
dist_genomes <- dist_sccg %>%
    left_join(dist_jaccard) %>%
    left_join(dist_ani) %>%
    left_join(dist_kmer)

nrow(dist_genomes) # 36^2 = 1296
write_csv(dist_genomes, paste0(folder_data, "genomics_analysis/distances/dist_genomes.csv"))


# Mantel test ----
convert_to_distm <- function (dist_genomes, d_i) {
    #' This function convert the long format distances into wide matrix format
    select(dist_genomes, genome_id1, genome_id2, {{d_i}}) %>%
    pivot_wider(names_from = genome_id2, values_from = {{d_i}}) %>%
        column_to_rownames(var = "genome_id1") %>%
        as.dist()
}
tidy_mantel <- function(mod) tibble(corr = mod$statistic, p = mod$signif)


list_d <- c("d_sccg", "d_jaccard", "d_ani", "d_kmer")
tb <- crossing(d1 = ordered(list_d, list_d), d2 = ordered(list_d, list_d)) %>%
    filter(d1 < d2) %>%
    rowwise() %>%
    mutate(
        mat1 = list(convert_to_distm(dist_genomes, d1)),
        mat2 = list(convert_to_distm(dist_genomes, d2)),
        mod_mantel = list(mantel(mat1, mat2, method = "pearson", permutations = 9999)),
        mod_tidied = list(tidy_mantel(mod_mantel))
    )

tb_mantel <- tb %>%
    unnest(mod_tidied) %>%
    select(d1, d2, corr, p)

write_csv(tb_mantel, paste0(folder_data, "genomics_analysis/distances/tb_mantel.csv"))
