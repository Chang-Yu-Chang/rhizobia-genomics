#' This script clean the gene content data

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

genomes <- read_csv(paste0(folder_data, "temp/00-genomes.csv")) 
#pa <- read_csv(paste0(folder_data, "genomics/pangenome/panaroo/gene_presence_absence.csv"))
pa <- read_delim(paste0(folder_data, "genomics/pangenome/panaroo/gene_presence_absence.Rtab"))
pa <- pa %>% clean_names()


# Compute pairwise distance in gene presence and absence 
pat <- t(as.matrix(pa[,-1]))

dim(pat)
rownames(pat) <- colnames(pa)[-1]
jdm <- proxy::dist(pat, method = "Jaccard")

dist_jaccard <- jdm %>%
    as.matrix() %>% as_tibble() %>%
    mutate(genome_id1 = colnames(.)) %>%
    pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "d_jaccard")

write_csv(dist_jaccard, paste0(folder_data, "temp/13-dist_jaccard.csv"))

# Compute genomic fluidity
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