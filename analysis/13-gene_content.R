#' This script clean the gene content data

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

genomes <- read_csv(paste0(folder_data, "temp/00-genomes.csv")) 
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
#pa <- read_csv(paste0(folder_data, "genomics/pangenome/panaroo/gene_presence_absence.csv"))
pa <- read_delim(paste0(folder_data, "genomics/pangenome/panaroo/gene_presence_absence.Rtab"))
pa <- pa %>% clean_names()

# Transpose the gene presence-absence table
gpa <- pa %>%
    pivot_longer(cols = -gene) %>%
    pivot_wider(names_from = gene, values_from = value, values_fill = 0) 

write_csv(gpa, paste0(folder_data, "temp/13-gpa.csv"))

# 1. Compute jaccard distance in gene presence and absence 
pat <- t(as.matrix(pa[,-1]))
dim(pat) # 41 x 31964
rownames(pat) <- colnames(pa)[-1]
jdm <- proxy::dist(pat, method = "Jaccard")

dist_jaccard <- jdm %>%
    as.matrix() %>% as_tibble() %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_jaccard")

write_csv(dist_jaccard, paste0(folder_data, "temp/13-dist_jaccard.csv"))

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
    pa_i <- pa[,c(dist_fluidity$genome_id1[i], dist_fluidity$genome_id2[i])] 
    pa_i <- pa_i[pa_i[,dist_fluidity$genome_id1[i]] == 1 | pa_i[,dist_fluidity$genome_id2[i]] == 1,]
    dist_fluidity$d_fluidity[i] <- fluidity(pa_i)
}

write_csv(dist_fluidity, paste0(folder_data, "temp/13-dist_fluidity.csv"))


# 3. Compute jaccard by contigs
# List of all genes in all genomes
pd <- read_csv(paste0(folder_data, "genomics/pangenome/panaroo/gene_data.csv"))
pd <- clean_names(pd) %>%
    select(genome_id = gff_file, contig_id = scaffold_name, annotation_id, description) %>%
    mutate(annotation_id = str_remove(annotation_id, ".+/genomes/"))

gpa <- read_csv(paste0(folder_data, "genomics/pangenome/panaroo/gene_presence_absence.csv")) 
gpa <- gpa %>%
    clean_names() %>%
    select(starts_with("g"), em1021, em1022, usda1106, wsm419) %>%
    pivot_longer(-gene, names_to = "genome_id", values_to = "annotation_id", values_drop_na = T) %>%
    mutate(annotation_id = str_remove(annotation_id, ".+/genomes/")) %>%
    select(genome_id, annotation_id, gene) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id, annotation_id)

pa_c <- gpa %>%
    left_join(pd) %>%
    #' The actual gene number in gpa is correct, but there are some genes in gpa that are not present in pd
    #' So I am not able to assign those genes to contigs
    drop_na(contig_id) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id)) 

nrow(pa_c) # 284528

# 
pa_ct <- pa_c %>%
    select(contig_id, gene) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = contig_id, values_from = value, values_fill = 0) # The number of total genes in the pangenome is correct
pa_ct <- t(as.matrix(pa_ct[,-1]))
dim(pa_ct) # 224 x 31964
jdm <- proxy::dist(pa_ct, method = "Jaccard")

dist_jaccard_contigs <- jdm %>%
    as.matrix() %>% as_tibble() %>%
    mutate(contig_id1 = colnames(.)) %>%
    pivot_longer(-contig_id1, names_to = "contig_id2", values_to = "d_jaccard")


write_csv(dist_jaccard_contigs, paste0(folder_data, "temp/13-dist_jaccard_contigs.csv"))
