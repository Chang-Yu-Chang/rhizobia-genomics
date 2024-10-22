#' This script generates the list of gene names for searching on Uniprot

library(tidyverse)
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

#set_name = "elev_med"
set_name = "urbn_mel"
tt <- read_gpas(set_name)

xx <- tt$gpar %>%
    mutate(across(matches("^g\\d"), ~str_remove_all(.x, "/Users/cychang/Dropbox/lab/rhizobia-genomics/data/genomics/fasta/genomes/"))) %>%
    pivot_longer(cols = -c(gene, non_unique_gene_name, annotation)) %>%
    separate_rows(value, sep = ";") %>%
    drop_na(value)

nrow(xx) # 73661

# Remove genes with no good annotation -> there are ~4647 unique genes have annotation among ~10k pangenes
yy <- xx %>%
    filter(!(str_detect(gene, "group") & is.na(non_unique_gene_name) & annotation == "hypothetical protein")) %>%
    distinct(gene, non_unique_gene_name, annotation, name, .keep_all = T) %>%
    distinct(gene, .keep_all = T)

list_unique_genes <- yy %>%
    filter(!str_detect(gene, "group")) %>%
    select(gene) %>%
    separate_rows("gene", sep = "~~~") %>%
    mutate(gene = str_remove(gene, "_\\d")) %>%
    unique() %>%
    arrange(gene)
nrow(list_unique_genes)
write_csv(list_unique_genes, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_unique_genes.csv"))

#list_unique_genes$gene %>% cat() # All 2000 genes with known names

#list_gene$gene %>% cat(sep = '", "') # All 2000 genes with known names
#GetProteinGOInfo("Q05572")

