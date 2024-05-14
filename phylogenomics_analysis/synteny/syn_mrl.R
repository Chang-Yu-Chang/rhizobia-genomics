#' This script follows Zhao 2021 pipeline for Syn-MRL. It includes four steps
#' 1. phylogenomic synteny network construction
#'  1a. all vs all reciprocal sequence similarity search for all annotated proteomes using DIAMOND
#'  1b. pairwise synteny block detection using MCScanX
#' 2. network clustering
#' 3. matrix representation
#' 4. ML based tree inference

renv::load()
library(tidyverse)
library(syntenet)
source(here::here("metadata.R"))

# Read the data ----
aa <- fasta2AAStringSetlist(paste0(folder_data, "genomics/faa/genomes"))
gr <- gff2GRangesList(paste0(folder_data, "genomics/gff/genomes"))
length(aa)
length(gr)

## Subset the list
aa <- aa[names(aa) %in% isolates$genome_id[-23]]
length(aa)
names(aa)
gr <- gr[names(gr) %in% isolates$genome_id[-23]]
length(gr)
names(gr)

## Clean the gene ids
clean_aa_names <- function (aastrings) {
    names(aastrings) <- names(aastrings) %>%
        str_remove(".+/genomes/") %>%
        str_remove(" .+")
    return(aastrings)
}
clean_gr_ids <- function (granges) {
    granges$ID <- granges$ID %>% str_remove(".+/genomes/")
    granges$Name <- granges$ID
    granges$gene_id <- granges$ID

    genome_id <- str_remove(granges$ID[1], "\\..+")
    new_seqnames <- paste0(genome_id, "_", seqnames(granges))
    seqlevels(gr) <- levels(new_seqnames)
    seqnames(granges) <- new_seqnames
    return(granges[granges$type == "CDS",c("source", "type", "score", "phase", "ID", "Name", "gene_id")])
}

aa <- lapply(aa, clean_aa_names)
gr <- lapply(gr, clean_gr_ids)

##

names(aa)
names(gr)
aa$g10
gr$g10
sapply(aa, length)
sapply(gr, length)

check_input(aa, gr)
process_input(aa, gr)
length(aa$g10)
length(gr$g10)



