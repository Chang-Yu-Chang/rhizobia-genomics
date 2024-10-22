#' This script performs the  the GO terms
#' On Uniprot ID mapping https://www.uniprot.org/id-mapping
#' The list of unique genes were searched against four species level taxa IDs:
#' tax_id 382 Rhizobium meliloti (species)
#' tax_id 110321 Sinorhizobium medicae (species)
#' tax_id 562  Escherichia coli (species)
#' tax_id 287 Pseudomonas aeruginosa (species)

library(tidyverse)
library(janitor)
#library(GO.db) # A set of annotation maps describing the entire Gene Ontology
library(topGO) # for GO term enrichment
library(ggsci)
library(cowplot)
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

#set_name = "elev_med"
set_name = "urbn_mel"
tt <- read_gpas(set_name)
ff <- read_fsts(set_name)
gg <- read_gos(set_name)


# The list of all genes
go_ids <- gg$gg_all %>%
    distinct(from, .keep_all = T) %>%
    mutate(go_id = str_extract_all(gene_ontology_go, "GO:\\d+"))
    #dplyr::select(from, entry, organism, go_id, gene_ontology_go) %>%
    #unnest(cols = go_id) %>%
    #mutate(go_id = )
    # dplyr::select(-gene_ontology_go) %>%
    # mutate(gop_id = map(go_id, ~ggp[which(names(ggp) == .x)]))

# Mapping object (a R list) required for topGO
gene2go <- tt$gpa %>%
    # Clean up gene names
    mutate(gene = str_split(gene, "~~~")) %>%
    unnest(gene) %>%
    mutate(from = str_remove(gene, "_\\d")) %>%
    filter(!str_detect(gene, "group")) %>%
    left_join(go_ids) %>%
    drop_na(gene_ontology_go) %>%
    dplyr::select(gene, go_id) %>%
    deframe()

# Genes with top 5% Fst
top_gene_fst <- ff$per_gene_fst %>%
    left_join(ff$gene_lengths) %>%
    # Remove unannotated genes
    filter(!str_detect(gene, "group")) %>%
    arrange(desc(fst)) %>%
    # Choose the top 5% genes
    mutate(tops = ifelse(fst > quantile(fst, 0.95), 1, 0)) %>%
    # Clean up gene names
    mutate(gene = str_split(gene, "~~~")) %>%
    unnest(gene)
    #mutate(from = str_remove(gene, "_\\d")) %>%


gene_interested <- top_gene_fst %>%
    mutate(tops = factor(tops, c(1, 0))) %>%
    dplyr::select(gene, tops) %>%
    deframe()

# Create the topGOdata object
oo <- c("BP", "CC", "MF")

list_godata <- list()
list_fisher <- list()
list_tables <- list()

for (o in oo) {
    godata <- new("topGOdata",
                  ontology = o,
                  allGenes = gene_interested,
                  nodeSize = 15,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2go
    )
    list_godata[[o]] <- godata
    list_fisher[[o]] <- runTest(godata, algorithm = "classic", statistic = "fisher")
    list_tables[[o]] <- GenTable(list_godata[[o]], classicFisher = list_fisher[[o]])

}

goenrich <- bind_rows(list_tables, .id = "Category") %>%
    clean_names %>%
    mutate(goidterm = paste0(go_id, " ", term))


p1 <- goenrich %>%
    mutate(goidterm = factor(goidterm, rev(goenrich$goidterm))) %>%
    ggplot() +
    geom_col(aes(x = goidterm, y = annotated, fill = category)) +
    coord_flip() +
    scale_fill_aaas() +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(.9,.9),
        legend.background = element_rect(color = "black"),
        axis.text.y.left = element_text(hjust = 0)
    ) +
    guides() +
    labs(x = "", title = set_name)

p2 <- goenrich %>%
    mutate(goidterm = factor(goidterm, rev(goenrich$goidterm))) %>%
    ggplot() +
    geom_col(aes(x = goidterm, y = -log(as.numeric(classic_fisher), 10), fill = category)) +
    geom_hline(yintercept = -log(0.05, 10)) +
    coord_flip() +
    scale_fill_aaas() +
    theme_bw() +
    theme(
        axis.text.y = element_blank(),
        legend.position = "none"
    ) +
    guides() +
    labs(x = "", y = "-log(p)")

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", rel_widths = c(1, .8))
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-05-go.png"), p, width = 15, height = 8)


# o = "BP"
# GenTable(list_godata[[o]], classicFisher = list_fisher[[o]])
# o = "CC"
# GenTable(list_godata[[o]], classicFisher = list_fisher[[o]])
# o = "MF"
# GenTable(list_godata[[o]], classicFisher = list_fisher[[o]])






















