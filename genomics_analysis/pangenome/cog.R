#' This script finds the GO terms for the top gene clusters

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- arrange(isolates, site_group)
read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/spa.csv"))

    return(list(gpa = gpa, gene_order = gene_order, gpatl = gpatl, gpacl = gpacl, gd = gd, sml = sml, list_sccg = list_sccg, spa = spa))
}

# Download NIH COG list
# https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab.txt

cog24 <- read_delim( "~/Downloads/cog-24.def.tab.txt", delim = "\t", col_names = c("COG ID", "COG functional category", "COG name", "Gene name", "Functional pathway", "PubMed ID", "PDB ID"))
cog24 <- cog24 %>% clean_names()

#
set_name = "urbn_mel"
#set_name = "elev_med"
top_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name, "/top_gene_fst.csv")) %>% filter(!str_detect(gene, "group"))
list_tg_fst <- top_gene_fst$gene %>% str_split("~~~") %>% unlist()
#cat(paste(gene_names, collapse = " "))
top_gene_or <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/top_gene_or.csv")) %>% filter(!str_detect(gene, "group"))
list_tg_or <- top_gene_or$gene %>% str_split("~~~") %>% unlist()


# Need the prokka output tsv table for COG id
#genome_id <- paste0("g", c(4,5,6,8,9,11,13,16,17,19))
genome_id <- paste0("g", c(21:27, 31:37, 39, 41, 43))
list_gff <- list()
for (gi in genome_id) {
    list_gff[[gi]] <- read_delim(paste0(folder_data, "genomics/annotation/", gi,"/annotated.tsv"), delim = "\t", skip = 0) %>%
        clean_names() %>%
        rename(cog_id = cog)
}

gffs <- bind_rows(list_gff, .id = "genome_id") #%>% distinct(gene, .keep_all = T)
gffs %>%
    filter(str_detect(gene, "abo")) %>%
    view

cog_tg_fst <- gffs %>%
    filter(gene %in% list_tg_fst) %>%
    left_join(cog24)

table(cog_tg_fst$cog_functional_category)

cog_tg_or <- gffs %>%
    filter(gene %in% list_tg_or) %>%
    left_join(cog24)

table(cog_tg_or$cog_functional_category)










