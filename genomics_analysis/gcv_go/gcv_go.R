#' This script performs the GO enrichment analysis

library(tidyverse)
library(janitor)
library(topGO) # for GO term enrichment
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
read_gcv_fsts <- function (set_name) {
    per_acce_fst <- read_csv(paste0(folder_data, "genomics_analysis/gcv_fst/", set_name,"/per_acce_fst.csv"))
    per_genome_fst <- read_csv(paste0(folder_data, "genomics_analysis/gcv_fst/", set_name,"/per_genome_fst.csv"))
    return(list(per_acce_fst = per_acce_fst, per_genome_fst = per_genome_fst))
}
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
read_gos_generic <- function (set_name) {
    generic <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_generic")) %>% clean_names()
    return(generic)
}
make_gene2go <- function (tt, go_ids) {
    tt$gpa %>%
        # Clean up gene names
        mutate(gene = str_split(gene, "~~~")) %>%
        unnest(gene) %>%
        mutate(from = str_remove(gene, "_\\d")) %>%
        filter(!str_detect(gene, "group")) %>%
        left_join(go_ids) %>%
        drop_na(gene_ontology_go) %>%
        dplyr::select(gene, go_id) %>% # Make sure the paralogs are used
        deframe()
}
make_top_genes_bygene <- function (ff, pp = 0.95) {
    #' pp = 0.95 means selecting the top 5%
    ff$per_acce_fst %>%
        # Remove unannotated genes
        filter(!str_detect(gene, "group")) %>%
        arrange(desc(Gprime_st)) %>%
        drop_na(Gprime_st) %>%
        # Choose the top 5% genes
        mutate(tops = ifelse(Gprime_st >= quantile(Gprime_st, pp), 1, 0)) %>%
        # Clean up gene names
        mutate(gene = str_split(gene, "~~~")) %>%
        unnest(gene) %>%
        distinct(gene, .keep_all = T)
}
make_ginterest <- function (top_genes) {
    top_genes %>%
        mutate(tops = factor(tops, c(1, 0))) %>%
        dplyr::select(gene, tops) %>%
        deframe()
}
go_wrapper <- function (set_name, ginterest, oo, suffix = "", gene2go) {
    #' A wrapper function for go enrichment analysis
    list_godata <- list()
    list_fisher <- list()
    list_tables <- list()
    for (o in oo) {
        # Create the topGOdata object for each category
        godata <- new(
            "topGOdata",
            ontology = o,
            allGenes = ginterest,
            nodeSize = 10,
            annot = annFUN.gene2GO,
            gene2GO = gene2go
        )
        list_godata[[o]] <- godata
        list_fisher[[o]] <- runTest(godata, algorithm = "classic", statistic = "fisher")
        list_tables[[o]] <- GenTable(list_godata[[o]], classicFisher = list_fisher[[o]])
    }
    goenrich <- bind_rows(list_tables, .id = "Category") %>% clean_names
    write_csv(goenrich, paste0(folder_data, "genomics_analysis/gcv_go/", set_name, "/goenrich", suffix,".csv"))
}
total_wrapper <- function (set_name) {
    tt <- read_gpas(set_name)
    ff <- read_gcv_fsts(set_name)
    gg <- read_gos(set_name)
    gg$generic <- read_gos_generic(set_name) # GO ids from genric search

    # The list of GO ids
    go_ids <- gg$gg_all %>%
        bind_rows(gg$generic) %>%
        drop_na(gene_ontology_go) %>%
        distinct(from, .keep_all = T) %>%
        mutate(go_id = str_extract_all(gene_ontology_go, "GO:\\d+"))

    round(nrow(go_ids)/nrow(gg$list_unique_genes), 3) # % of genes that have GO ids

    # Mapping object (a R list) required for topGO. It the gene universe
    gene2go <- make_gene2go(tt, go_ids)
    length(gene2go) # number of clusters that have GO terms in the universe

    # Genes with top 5% Fst
    top_genes_bygene <- make_top_genes_bygene(ff) # by gene wide Fst
    write_csv(top_genes_bygene, paste0(folder_data, "genomics_analysis/gcv_go/", set_name, "/top_genes_bygene.csv"))

    # Make gene of interest
    ginterest_bygene <- make_ginterest(top_genes_bygene)
    length(ginterest_bygene)
    table(ginterest_bygene)

    oo <- c("BP", "CC", "MF")
    go_wrapper(set_name, ginterest_bygene, oo, suffix = "_bygene", gene2go)
    #go_wrapper(set_name, ginterest_bysnp, oo, suffix = "_bysnp")
}

total_wrapper("elev_med")
total_wrapper("urbn_mel")
