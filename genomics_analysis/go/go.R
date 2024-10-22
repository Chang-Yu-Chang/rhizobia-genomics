#' This script performs the  the GO terms
#' On Uniprot ID mapping https://www.uniprot.org/id-mapping
#'
#' The list of unique genes were searched against four species level taxa IDs:
#' tax_id 382 Rhizobium meliloti (species)
#' tax_id 110321 Sinorhizobium medicae (species)
#' tax_id 562  Escherichia coli (species)
#' tax_id 287 Pseudomonas aeruginosa (species)

library(tidyverse)
library(janitor)
#library(GO.db) # A set of annotation maps describing the entire Gene Ontology
library(topGO) # for GO term enrichment
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpar <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpar.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/spa.csv"))

    return(list(gpa = gpa, gene_order = gene_order, gpar = gpar, gpatl = gpatl, gpacl = gpacl, gd = gd, sml = sml, list_sccg = list_sccg, spa = spa))
}
read_fsts <- function (set_name) {
    per_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_gene_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))
    return(list(per_gene_fst = per_gene_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths))
}
read_gos <- function (set_name) {
    # List of unique genes used for Uniprot ID mapping
    list_unique_genes <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_unique_genes.csv"))
    # tax_id 382 Rhizobium meliloti (species)
    to_med <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_med.tsv")) %>% clean_names()
    # tax_id 110321 Sinorhizobium medicae (species)
    to_mel <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_mel.tsv")) %>% clean_names()
    # tax_id 562  Escherichia coli (species)
    to_eco <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_eco.tsv")) %>% clean_names()
    # tax_id 287 Pseudomonas aeruginosa (species)
    to_par <- read_tsv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/", set_name, "_to_par.tsv")) %>% clean_names()
    # Drop the entries with no GO description
    gg_all <- bind_rows(to_med, to_mel, to_par, to_eco) %>% drop_na(gene_ontology_go)

    return(list(list_unique_genes = list_unique_genes, to_med = to_med, to_mel = to_mel, to_eco = to_eco, to_par = to_par, gg_all = gg_all))
}

make_gene2go <- function (gpa, go_ids) {
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
    ff$per_gene_fst %>%
        left_join(ff$gene_lengths) %>%
        # Remove unannotated genes
        filter(!str_detect(gene, "group")) %>%
        arrange(desc(fst)) %>%
        # Choose the top 5% genes
        mutate(tops = ifelse(fst > quantile(fst, pp), 1, 0)) %>%
        # Clean up gene names
        mutate(gene = str_split(gene, "~~~")) %>%
        unnest(gene) %>%
        distinct(gene, .keep_all = T)
}
make_top_genes_bysnp <- function (ff, pp = 0.95) {
    ff$per_locus_fst %>%
        # Remove unannotated genes
        filter(!str_detect(gene, "group")) %>%
        left_join(ff$gene_lengths) %>%
        arrange(desc(fst)) %>%
        mutate(tops = ifelse(fst > quantile(fst, pp), 1, 0)) %>%
        # Clean up gene names
        mutate(gene = str_split(gene, "~~~")) %>%
        unnest(gene) %>%
        arrange(desc(tops), gene) %>%
        distinct(gene, tops, .keep_all = T)
}
make_top_genes_bygene_byreplicon <- function (ff, tt, pp = 0.95) {
    #' pp = 0.95 means selecting the top 5%
    # Use the first genome
    gene_replicon <- filter(tt$gpacl, genome_id == tt$gpacl$genome_id[1]) %>%
        filter(!str_detect(gene, "group")) %>%
        dplyr::select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))

    top_genes <- ff$per_gene_fst %>%
        left_join(ff$gene_lengths) %>%
        # Remove unannotated genes
        filter(!str_detect(gene, "group")) %>%
        # Compute by replicon
        left_join(gene_replicon) %>%
        filter(!is.na(replicon_type), replicon_type != "others") %>%
        group_by(replicon_type) %>%
        arrange(replicon_type, desc(fst)) %>%
        # Choose the top 5% genes
        mutate(tops = ifelse(fst > quantile(fst, pp), 1, 0)) %>%
        # Clean up gene names
        mutate(gene = str_split(gene, "~~~")) %>%
        unnest(gene)

    return(top_genes)
    # check
    # top_genes %>%
    #     group_by(replicon_type, tops) %>%
    #     count()
}
make_top_genes_bygene_bysnp <- function (ff, tt, pp = 0.95) {
    #' pp = 0.95 means selecting the top 5%
    # Use the first genome
    gene_replicon <- filter(tt$gpacl, genome_id == tt$gpacl$genome_id[1]) %>%
        filter(!str_detect(gene, "group")) %>%
        dplyr::select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))


    ff$per_locus_fst %>%
        left_join(ff$gene_lengths) %>%
        # Remove unannotated genes
        filter(!str_detect(gene, "group")) %>%
        # Compute by replicon
        left_join(gene_replicon) %>%
        filter(!is.na(replicon_type), replicon_type != "others") %>%
        #group_by(replicon_type) %>%
        arrange(replicon_type, desc(fst)) %>%
        # Choose the top 5% genes
        mutate(tops = ifelse(fst > quantile(fst, pp), 1, 0)) %>%
        # Clean up gene names
        mutate(gene = str_split(gene, "~~~")) %>%
        unnest(gene) %>%
        arrange(replicon_type, desc(tops), gene) %>%
        distinct(gene, tops, .keep_all = T)


    return(top_genes)
    # check
    top_genes %>%
        group_by(replicon_type, tops) %>%
        count()
}
make_ginterest <- function (top_genes) {
    top_genes %>%
        mutate(tops = factor(tops, c(1, 0))) %>%
        dplyr::select(gene, tops) %>%
        deframe()
}
go_wrapper <- function (ginterest, oo, suffix = "") {
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
    write_csv(goenrich, paste0(folder_data, "genomics_analysis/go/", set_name, "/goenrich", suffix,".csv"))
}


set_name = "elev_med"
#set_name = "urbn_mel"
tt <- read_gpas(set_name)
ff <- read_fsts(set_name)
gg <- read_gos(set_name)

# The list of all genes in this gradient pangenome
go_ids <- gg$gg_all %>%
    distinct(from, .keep_all = T) %>%
    mutate(go_id = str_extract_all(gene_ontology_go, "GO:\\d+"))

# Mapping object (a R list) required for topGO. It the gene universe
gene2go <- make_gene2go(tt$gpa, go_ids)
length(gene2go) # number of all genes in the universe

# Genes with top 5% Fst
top_genes_bygene <- make_top_genes_bygene(ff) # by gene wide Fst
top_genes_bysnp <- make_top_genes_bysnp(ff) # by snp Fst
#top_genes_bygene_byreplicon <- make_top_genes_bygene_byreplicon(ff, tt) # by gene by replicon

# Make gene of interest
#gene_interested <- make_gene_interested(top_gene_fst)
ginterest_bygene <- make_ginterest(top_genes_bygene)
length(ginterest_bygene)
table(ginterest_bygene)
ginterest_bysnp <- make_ginterest(top_genes_bysnp)
length(ginterest_bysnp)
table(ginterest_bysnp)

oo <- c("BP", "CC", "MF")
go_wrapper(ginterest_bygene, oo, suffix = "_bygene")
go_wrapper(ginterest_bysnp, oo, suffix = "_bysnp")



