#' This script cleans the tables of genetic distance and join them

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("metadata.R"))

folder_genomes <- paste0(folder_genomics, "fasta/genomes/")
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))

# 1. ANI
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
    filter(genome_id1 <= genome_id2) %>%
    select(genome_id1, genome_id2, d_ani)
#write_csv(dist_ani, paste0(folder_data, "temp/19-dist_ani.csv"))

# 2. k-mers
dist_kmer <- read_delim(paste0(folder_data, "genomics/kmer/genomes/kmer.txt"))
list_sigs <- read_delim(paste0(folder_data, "genomics/kmer/genomes/list_sigs.txt"), delim = "\t", col_names = F)

dist_kmer <- dist_kmer %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_kmer")  %>%
    mutate(genome_id1 = str_remove(genome_id1, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id2 = str_remove(genome_id2, ".+/genomes/") %>% str_remove(".fasta")) %>%
    mutate(genome_id1 = ordered(genome_id1, genomes$genome_id), genome_id2 = ordered(genome_id2, genomes$genome_id)) %>%
    arrange(genome_id1, genome_id2) %>%
    filter(genome_id1 <= genome_id2)

#write_csv(dist_kmer, paste0(folder_data, "temp/19-dist_kmer.csv"))

# 3. gene content
## jaccard distance in gene presence and absence
pat <- as.matrix(gpat[,-1])
dim(pat) # 41 x 31964
rownames(pat) <- gpat$genome_id
jdm <- proxy::dist(pat, method = "Jaccard")

dist_jaccard <- jdm %>%
    as.matrix() %>% as_tibble() %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_jaccard")

#write_csv(dist_jaccard, paste0(folder_data, "genomics_analysis/gene_content/dist_jaccard.csv"))

# 2. Compute genomic fluidity in pairs
fluidity <- function (tb) {
    uk = sum(tb[,1] == 1 & tb[,2] == 0)
    ul = sum(tb[,1] == 0 & tb[,2] == 1)
    mk = sum(tb[,1])
    ml = sum(tb[,2])
    gf <- (uk + ul) / (mk + ml)
    return(gf)
}

dist_fluidity <- dist_jaccard %>%
    select(genome_id1, genome_id2) %>%
    mutate(d_fluidity = NA)

for (i in 1:nrow(dist_fluidity)) {
    pa_i <- gpa[,c(dist_fluidity$genome_id1[i], dist_fluidity$genome_id2[i])]
    pa_i <- pa_i[pa_i[,dist_fluidity$genome_id1[i]] == 1 | pa_i[,dist_fluidity$genome_id2[i]] == 1,]
    dist_fluidity$d_fluidity[i] <- fluidity(pa_i)
}


# Join the distance tables
dist_genetics <- dist_ani %>%
    left_join(dist_kmer) %>%
    left_join(dist_jaccard) %>%
    left_join(dist_fluidity) %>%
    filter(genome_id1 != genome_id2)

nrow(dist_genetics) # choose(32,2)=496

write_csv(dist_genetics, paste0(folder_data, "genomics_analysis/dist_genetics.csv"))



# # 4. 825 core genes
# dist_scccg <- read_csv(paste0(folder_data, "temp/15-dist_sccg.csv"))
#dist_jaccard_contigs <- read_csv(paste0(folder_data, "temp/13-dist_jaccard_contigs.csv"))
# # 2.2 contig kmers
# dist_kmer_contigs <- read_delim(paste0(folder_data, "genomics/popgen/kmer_contigs/kmer.txt"))
# list_sigs <- read_delim(paste0(folder_data, "genomics/popgen/kmer_contigs/list_sigs.txt"), delim = "\t", col_names = F)
# list_sigs <- str_extract(list_sigs$X1, "g\\d+_contig_\\d+")
# length(list_sigs) # 228 contigs
# dist_kmer_contigs <- dist_kmer_contigs %>%
#     mutate(contig_id1 = colnames(.)) %>%
#     pivot_longer(-contig_id1, names_to = "contig_id2", values_to = "d_kmer")  %>%
#     mutate(contig_id1 = str_remove(contig_id1, ".+/contigs/") %>% str_remove(".fasta")) %>%
#     mutate(contig_id2 = str_remove(contig_id2, ".+/contigs/") %>% str_remove(".fasta")) %>%
#     mutate(genome_id1 = str_remove(contig_id1, "_contig_\\d"), genome_id2 = str_remove(contig_id2, "_contig_\\d")) %>%
#     mutate(contig_id1 = ordered(contig_id1, list_sigs), contig_id2 = ordered(contig_id2, list_sigs)) %>%
#     arrange(contig_id1, contig_id2) %>%
#     filter(contig_id1 <= contig_id2) %>%
#     select(starts_with("genome_id"), starts_with("contig"), d_kmer)
#
# nrow(dist_kmer_contigs) # choose(228,2) + 228 = 26106
#
#
#
# dist_genetics_contigs <- dist_kmer_contigs %>%
#     left_join(dist_jaccard_contigs)
#
# write_csv(dist_genetics_contigs, paste0(folder_data, "temp/19-dist_genetics_contigs.csv"))
#
#
#
#
#
#
