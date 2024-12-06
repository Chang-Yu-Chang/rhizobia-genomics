#' This script checks the gene names and annotations and output a list to be searched in uniprot

library(tidyverse)
library(janitor)
library(vegan)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

#set_name = "elev_med"
set_name = "urbn_mel"
tt <- read_gpas(set_name)

# Check the number of gene symbols -----

# How many clusters in the pangnome?
nt <- nrow(tt$gpa)
cat("\nTotal genes in pangenome:", nt, "/", nt, "(", round(nt/nt*100, 2),"%)")

# How many core genes/clusters?
nn <- tt$gpar %>% filter(no_isolates == max(no_isolates)) %>% nrow()
cat("\nCore genes in pangenome:", nn,  "/", nt,"(", round(nn/nt*100, 2),"%)")
# Single copy core genes
nn <- tt$gpar %>% filter(no_isolates == max(no_isolates)) %>% filter(avg_sequences_per_isolate == 1) %>% nrow
cat("\nSingle copy core genes in pangenome:", nn, "/", nt, "(", round(nn/nt*100, 2),"%)")

# How many gene clusters are annotated?
#' It's annotated when it meets at least one of the three criteria
#' 1. `gene` is not "group_XX"
#' 2. `non_unique_gene_name` is not NA
#' 3. `annotation` is not "hypothetical protein" nor "putative protein"

gparsub <- tt$gpar %>%
    filter(!str_detect(gene, "group") | !is.na(non_unique_gene_name) | !(annotation %in% c("hypothetical protein", "putative protein")))

nn <- nrow(gparsub) # number of annotated clusters
cat('\nAnnotation criteria\n1. `gene` is not "group_XX" \n2. `non_unique_gene_name` is not NA\n3. `annotation` is not "hypothetical protein" nor "putative protein"')
cat("\nAnnotated genes in pangenome:", nn, "/", nt, "(", round(nn/nt*100, 2),"%)")

# Among annotated genes, how many of them have a gene name/symbols?
#' The gene name is either in `gene` or `non_unique_gene_name`
tb1 <- gparsub %>%
    filter(!str_detect(gene, "group"))

tb2 <- gparsub %>%
    filter(str_detect(gene, "group")) %>%
    # The gene names is in non_unique_gene_name
    filter(!is.na(non_unique_gene_name)) %>%
    # Clean non_unique_gene_name
    mutate(
        non_unique_gene_name = non_unique_gene_name %>%
            str_remove("^;") %>%
            str_remove(";.+") %>%
            str_remove(";")
    )

nn1 <- nrow(tb1) # number of clusters with a gene name in column `gene`
nn2 <- nrow(tb2) # number of clusters with a gene name in column `non_unique_gene_name`
tb <- bind_rows(tb1, tb2); nn <- nrow(tb)
nrow(tb) / nrow(tt$gpar) # % of clusters with gene names
cat("\nGenes with gene symbols in pangenome:", nn, "/", nt, "(", round(nn/nt*100, 2),"%)")
cat("\nGenes with one gene symbol in column `gene`:", nn1, "/", nn, "(", round(nn1/nn*100, 2),"%)")
cat("\nGenes with one gene symbol in column `non_unique_gene_name`:", nn2, "/", nn, "(", round(nn2/nn*100, 2),"%)")

# Among these genes with names, how many of these genes have more than two gene names?
#' An example would be like hdfR_3~~~hdfR_2 or  ugpC_3~~~ugpC_13~~~btuD_3
#' This is a gene cluster, but each copy was annotated in different names in PROKKA
tb_tg <- tb %>% filter(str_detect(gene, "~~~"))
nn <- nrow(tb_tg)
cat("\nGenes with multiple gene symbols:", nn, "/", nt, "(", round(nn/nt*100, 2),"%)")
cat("\nAn example would be like hdfR_3~~~hdfR_2 or ugpC_3~~~ugpC_13~~~btuD_3")

# Among these genes with two or more gene names, how many have different gene names?
#' For instance, hdfR_3~~~hdfR_2 has the same gene name hdfR
#' whereas ugpC_3~~~ugpC_13~~~btuD_3 has two gene names ugpC and btuD

tb_tgd <- tb_tg %>%
    mutate(gene = str_remove_all(gene,"_\\d")) %>%
    mutate(gg = str_split(gene, "~~~")) %>%
    mutate(is_same_gene = map_lgl(gg, ~length(unique(.x)) ==1)) %>%
    filter(!is_same_gene)
nn <- nrow(tb_tgd)
cat("\nGenes with multiple gene symbols and these gene symbols are not all equal:", nn, "/", nt, "(", round(nn/nt*100, 2),"%)")
cat("\nFor example\n\thdfR_3~~~hdfR_2 has the same gene name hdfR\n\twhereas ugpC_3~~~ugpC_13~~~btuD_3 has two different gene symbols ugpC and btuD")

# Among these genes with two or more different gene names, how many have a concensus?
#' For instance, grsB~~~grsB~~~grsB~~~grsB~~~lgrD has a consensus on grsB
#' 	aaeA~~~yiaV is inconclusive
tb_drop <- tb_tgd %>%
    mutate(is_concluded = map_lgl(gg, ~length(unique(table(.x))) > 1)) %>%
    filter(!is_concluded)

nn <- nrow(tb_drop)
cat("\nGenes with multiple gene symbols, which are not all equal, and are not conclusive:", nn, "/", nt, "(", round(nn/nt*100, 2),"%)")
cat("\nFor instance\n\tgrsB~~~grsB~~~grsB~~~grsB~~~lgrD has a consensus on grsB\n\taaeA~~~yiaV is inconclusive")


# Write a csv of a list of unique gene names that will be used to find GO terms
#' This will be joining two sets of gene
#' 1. clusters with one gene name, and the names are in column `gene`
#' 2. clusters with one gene name, and the names are in column `non_unique_gene_name`
#' 3. clusters with two or more gene names, but with consensus
t1 <- tb %>%
    filter(!str_detect(gene, "~~~")) %>%
    filter(!str_detect(gene, "group")) %>%
    mutate(from = str_remove_all(gene,"_\\d+")) %>%
    select(from, everything())

t2 <- tb %>%
    filter(!str_detect(gene, "~~~")) %>%
    filter(str_detect(gene, "group")) %>%
    mutate(from = str_remove_all(non_unique_gene_name, "_\\d+|'")) %>%
    select(from, everything())

t3 <- tb_tgd %>%
    mutate(is_concluded = map_lgl(gg, ~length(unique(table(.x))) > 1)) %>%
    filter(is_concluded) %>%
    mutate(from = map_chr(gg, ~names(sort(table(.x), decreasing = T))[1])) %>%
    select(from, everything())

nn1 <- nrow(t1) # 1. clusters with one gene name, and the names are in column `gene`
nn2 <- nrow(t2) # 2. clusters with one gene name, and the names are in column `non_unique_gene_name`
nn3 <- nrow(t3) # 3. clusters with two or more gene names, but with consensus
tts <- bind_rows(t1, t2, t3) %>%
    distinct(from, .keep_all = T) %>%
    select(from, gene, non_unique_gene_name, annotation)
nns1 <- nn1+nn2+nn3
nns2 <- nrow(tts)

cat("\nGenes with known gene symbols: ", nns1, "/", nt, "(", round(nns1/nt*100, 2),"%)")
cat("\nAmmong these")
cat("\n\tClusters with one gene name, and the names are in column `gene`: ", nn1, "/", nns1, "(", round(nn1/nns1*100, 2),"%)")
cat("\n\tClusters with one gene name, and the names are in column `non_unique_gene_name`: ", nn2, "/", nns1, "(", round(nn2/nns1*100, 2),"%)")
cat("\n\tClusters with two or more gene names, but with consensus: ", nn3, "/", nns1, "(", round(nn3/nns1*100, 2),"%)")

cat("\n Number of unique gene symbols:", nns2, "from", nns1)
write_csv(tts, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_unique_genes.csv"))
paste(tts$from, collapse = " ") # paste this to ID mapping in Uniprot as shown below
#' On Uniprot ID mapping https://www.uniprot.org/id-mapping
#'
#' The list of unique genes were searched against four species level taxa IDs:
#' tax_id 382 Rhizobium meliloti (species)
#' tax_id 110321 Sinorhizobium medicae (species)
#' tax_id 380 Rhizobium fredii (species)
#' tax_id 562  Escherichia coli (species)
#' tax_id 287 Pseudomonas aeruginosa (species)
#' tax id 4932 Saccharomyces cerevisiae

# (DO THIS section after uniprot searches against individual taxa) Do a genric search on the genes that fail to find GO ----

read_gos <- function (set_name) {
    # List of unique genes used for Uniprot ID mapping
    list_unique_genes <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_unique_genes.csv"))
    # tax_id 382 Rhizobium meliloti (species)
    to_med <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_med.tsv")) %>% clean_names()
    # tax_id 110321 Sinorhizobium medicae (species)
    to_mel <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_mel.tsv")) %>% clean_names()
    # tax_id 380 Rhizobium fredii (species)
    to_fre <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_fre.tsv")) %>% clean_names()
    # tax_id 562  Escherichia coli (species)
    to_eco <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_eco.tsv")) %>% clean_names()
    # tax_id 287 Pseudomonas aeruginosa (species)
    to_par <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_par.tsv")) %>% clean_names()
    # Drop the entries with no GO description
    gg_all <- bind_rows(to_med, to_mel, to_fre, to_par, to_eco) %>% drop_na(gene_ontology_go)

    return(list(list_unique_genes = list_unique_genes, to_med = to_med, to_mel = to_mel, to_fre = to_fre, to_eco = to_eco, to_par = to_par, gg_all = gg_all))
}
gg <- read_gos(set_name)

# The list of GO ids for those gene symbols found in the five taxa
go_ids <- gg$gg_all %>%
    distinct(from, .keep_all = T) %>%
    mutate(go_id = str_extract_all(gene_ontology_go, "GO:\\d+"))

nrow(go_ids)/nrow(gg$list_unique_genes) # % of genes that have GO ids

list_for_generic <- gg$list_unique_genes %>%
    left_join(go_ids) %>%
    filter(is.na(entry)) %>%
    filter(!str_detect(annotation, "hypothetical")) %>%
    select(from, gene, non_unique_gene_name, annotation)

nrow(list_for_generic)
write_csv(list_for_generic, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_for_generic.csv"))

paste0(list_for_generic$from, collapse = " ")






