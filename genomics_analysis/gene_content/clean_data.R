#' This script cleans the gene presence/absence from panaroo output

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

# Gene presence-absence table -----
gpa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.Rtab")) %>%
    clean_names() # Rows are genes, columns are genomes
dim(gpa) # 26886 genes in union x 32-1 genomes
write_csv(gpa, paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))

# Transpose the gene presence-absence table
gpat <- gpa %>%
    pivot_longer(cols = -gene, names_to = "genome_id") %>%
    pivot_wider(names_from = gene, values_from = value, values_fill = 0)
dim(gpat) # 31 genomes x 26887 genes
write_csv(gpat, paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))

# Compute gene counts -----
gpatl <- gpat %>%
    pivot_longer(cols = -genome_id, names_to = "gene") %>%
    filter(value == 1)
gene_order <- gpatl %>%
    group_by(gene) %>%
    dplyr::count() %>%
    arrange(desc(n)) %>%
    pull(gene)
gpatl <- gpatl %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    mutate(gene = factor(gene, gene_order))
dim(gpatl) # 213135 3
length(gene_order) # 26886 genes order by copy number
write_csv(gpatl, paste0(folder_data, "genomics_analysis/gene_content/gpatl.csv"))
write_csv(tibble(gene = gene_order), paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv"))

# Structural variation -----
spa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/struct_presence_absence.Rtab")) %>%
    clean_names()
dim(spa) # 9775 structural variants x 32-1 genomes
write_csv(spa, paste0(folder_data, "genomics_analysis/gene_content/spa.csv"))


# Full directory for each gene -----
gpaf <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.csv")) %>%
    clean_names() %>%
    mutate(across(starts_with("g"), function (x) str_remove_all(x, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")))
dim(gpaf) # 26886 genes x 34 genomes. 1) gene name 2) non unique gene name 3) annotation
write_csv(gpaf, paste0(folder_data, "genomics_analysis/gene_content/gpaf.csv"))

# Genes and the contigs they are from
gd <- read_csv(paste0(folder_data, "genomics/pangenome/isolates/gene_data.csv")) %>%
    clean_names() %>%
    mutate(annotation_id = str_remove_all(annotation_id, "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/")) %>%
    rename(genome_id = gff_file, contig_id = scaffold_name) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id))
dim(gd) # 225019 (gene x genome x contig) x 8 rows
write_csv(gd, paste0(folder_data, "genomics_analysis/gene_content/gd.csv"))

# Genes with contig origin ----
## Clean the gpaf names
gpaf <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpaf.csv"))
dim(gpaf) # 26886 x 34
gpafl <- gpaf %>%
    select(-non_unique_gene_name) %>%
    pivot_longer(-c(gene, annotation), names_to = "genome_id", values_to = "annotation_id", values_drop_na = T) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id) %>%
    mutate(annotation_id = str_remove(annotation_id, "_stop"))
dim(gpafl) # 213135

## Unpack genes with multiple copies in one genome
split_colons <- function(string) unlist(str_split(string, ";"))
split_tildes <- function(string) unlist(str_split(string, "~~~"))
gpafl_multi <- gpafl %>%
    filter(str_detect(annotation_id, ";")) %>%
    rowwise() %>%
    mutate(annotation_id = list(split_colons(annotation_id))) %>%
    unnest(annotation_id)
gpafl <- gpafl %>%
    filter(!str_detect(annotation_id, ";")) %>%
    bind_rows(gpafl_multi)
dim(gpafl) # 221,882

# Clean the gd names
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv")) %>% select(contig_id, replicon_type) %>% drop_na
gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gd.csv"))
dim(gd) # 225019 (gene x genome x contig) x 8 rows
gd <- gd %>% select(genome_id, contig_id, annotation_id)

# Join the two tables
gpacl <- gpafl %>%
    arrange(genome_id, annotation_id) %>%
    left_join(gd) %>%
    # Remove duplicate of multi-copy genes on one contig
    distinct(gene, contig_id) %>%
    # Filter for chromosome and two plasmids
    left_join(contigs)

length(unique(gpacl$contig_id))

write_csv(gpacl, paste0(folder_data, "genomics_analysis/gene_content/gpacl.csv"))

if (F) {

# Filter for single copy genes ----
gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gd.csv")) %>% rename(gene = gene_name)
gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpacl.csv"))

tb <- gpacl %>%
    filter(str_detect(gene, "nifH")) %>%
    drop_na(replicon_type)
tb_multi <- tb %>%
    filter(str_detect(gene, "~~~")) %>%
    rowwise() %>%
    mutate(gene = list(split_tildes(gene))) %>%
    unnest(gene)
tb %>%
    filter(!str_detect(gene, "~~~")) %>%
    bind_rows(tb_multi) %>%
    #arrange(gene, contig_id) %>%
    left_join(gd) %>%
    drop_na(genome_id) %>%
    arrange(gene) %>%
    view




    left_join(gd)
    filter(str_detect(gene_name, "nif|fix|nod")) %>%
    mutate(gene_name_group = str_sub(gene_name, 1, 4)) %>%
    filter(genome_id == "g10", gene_name_group == "fixK") %>%
    #group_by(genome_id, gene_name_group) %>%
    #count
    #filter(str_detect(gene_name, "nifH")) %>%
    mutate(tt = str_count(dna_sequence)) %>%
    #filter(str_count(prot_sequence) > 100) %>%
    #mutate(gene_name = str_remove(gene_name, "_\\d+")) %>%
    #distinct(genome_id, contig_id, gene_name) %>%
    view






}
