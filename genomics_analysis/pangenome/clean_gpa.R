#' This script cleans the gene presence/absence from panaroo outputs
#' It works on three sets
#'  - 1. all 36 genomes including non-symbiotic ensifers
#'  - 2. 10 elevation medicae
#'  - 3. 17 urbanization meliloti
#'
#' For each set, there will be 6 csv files
#'  gpa, spa, gpatl, gene_order, gd, gpacl
#'

library(tidyverse)
library(janitor)
library(vegan)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

clean_gpa <- function (gpa_file) {
    #' This function cleans the raw gpa file
    # Gene presence-absence table
    gpa <- read_delim(gpa_file, show_col_types = F) %>% clean_names() # Rows are genes, columns are genomes
    #cat(dim(gpa)) # 26504 genes in union x 37-1 genomes
    return(gpa)
}
get_sccg <- function (gpa_csv) {
    #gpa_csv = paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/gpa.csv")
    xx <- read_csv(gpa_csv) %>%
        clean_names()

    yy <- xx[,c(1, 4:ncol(xx))]
    zz <- yy %>%
        pivot_longer(cols = -gene, values_drop_na = T) %>%
        mutate(value = str_replace_all(value, "/Users/cychang/Dropbox/lab/rhizobia-genomics/data/genomics/fasta/genomes/", ""))
        #mutate(is_paralog = str_detect(value, ";") | str_detect(value, "refound"))
    aa <- zz %>%
        filter(!str_detect(value, "refound"), !str_detect(value, ";")) %>%
        pivot_wider(id_cols = gene, names_from = name, values_from = value) %>%
        arrange(gene) %>%
        # Remove accessory genes
        drop_na()

    # Remove genes with paralogs
    #list_sccg <- aa$gene[apply(aa[,-1], 1, sum) == 0]
    list_sccg <- aa$gene
    length(list_sccg)
    return(tibble(gene = list_sccg))
}
make_long_gpa <- function (gpa) {
    # Transpose the gene presence-absence table
    gpat <- gpa %>%
        pivot_longer(cols = -gene, names_to = "genome_id") %>%
        pivot_wider(names_from = gene, values_from = value, values_fill = 0)
    #dim(gpat) # 36 genomes x 26505 genes
    gpatl <- gpat %>%
        pivot_longer(cols = -genome_id, names_to = "gene") %>%
        filter(value == 1)
    gene_order <- gpatl %>%
        group_by(gene) %>%
        dplyr::count() %>%
        arrange(desc(n))
    #dim(gpatl) # 244960 3
    return(list(gpatl = gpatl, gene_order = gene_order))
}
clean_spa <- function (spa_file) {
    spa <- read_delim(spa_file) %>% clean_names()
    return(spa)
}
clean_gd <- function (gd_file) {
    gd <- read_csv(gd_file) %>%
        clean_names() %>%
        #mutate(annotation_id = str_replace(annotation_id, ".*/g","g")) %>%
        rename(genome_id = gff_file, contig_id = scaffold_name, gene = gene_name) %>%
        mutate(contig_id = paste0(genome_id, "_", contig_id)) %>%
    return(gd)
}
clean_gpaf <- function (gpaf_file) {
    gpaf <- read_csv(gpaf_file) %>%
        clean_names()
    #mutate(across(starts_with("g"), function (x) x %>% str_split(";") %>% unlist %>% str_replace(".*/g", "g") %>% str_c(collapse = ";")))
    #mutate(across(starts_with("g"), function (x) str_replace(x, ".*/g", "g")))
    return(gpaf)
}
make_long_gpaf <- function (gpaf) {
    gpafl <- gpaf %>%
        select(-non_unique_gene_name) %>%
        pivot_longer(-c(gene, annotation), names_to = "genome_id", values_to = "annotation_id", values_drop_na = T) %>%
        mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
        arrange(genome_id) %>%
        mutate(annotation_id = str_remove(annotation_id, "_stop"))

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
    dim(gpafl) # 255485
    return(gpafl)
}
make_long_gpac <- function (gpafl, contigs, gd) {
    # contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv")) %>% select(contig_id, replicon_type) %>% drop_na
    # gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gd.csv"))
    gd <- gd %>% select(genome_id, contig_id, annotation_id)
    gpacl <- gpafl %>%
        arrange(genome_id, annotation_id) %>%
        left_join(gd) %>%
        # Remove duplicate of multi-copy genes on one contig
        distinct(gene, contig_id, .keep_all = T) %>%
        # Filter for chromosome and two plasmids
        left_join(contigs)
    return(gpacl)
}
compute_bc_dist <- function (gpa) {
    gene_data_matrix <- t(as.matrix(gpa[,-1]))
    bcs <- vegdist(gene_data_matrix, method = "bray")
    similarity_matrix <- 1 - as.matrix(bcs)
    sml <- similarity_matrix %>%
        as_tibble() %>%
        mutate(genome_id1 = names(gpa)[-1]) %>%
        pivot_longer(cols = -genome_id1, names_to = "genome_id2", values_to = "bray_curtis_similarity")
    return(sml)
}



clean_all <- function (set_name, contigs) {
    #' A wrapper function performing all functions above
    dir_path <- paste0(folder_data, "genomics_analysis/gene_content/", set_name)
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE) # Create the directory (and parent directories, if needed)

    # Gene presence-absence table. Rows are genes, columns are genomes
    gpa <- clean_gpa(paste0(folder_data, "genomics/pangenome/", set_name,"/gene_presence_absence.Rtab"))
    dim(gpa) # 26504 genes in union x 37-1 genomes
    write_csv(gpa, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/gpa.csv"))

    # Get the list of single-copy core genes
    list_sccg <- get_sccg(paste0(folder_data, "genomics/pangenome/", set_name,"/gene_presence_absence.csv"))
    nrow(list_sccg)
    write_csv(list_sccg, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/list_sccg.csv"), col_names = F)

    # Compute the pairwise bray-curtis simularity based on gene presence absence
    sml <- compute_bc_dist(gpa)
    write_csv(sml, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/sml.csv"))

    # Structural variation
    spa <- clean_spa(paste0(folder_data, "genomics/pangenome/", set_name,"/struct_presence_absence.Rtab"))
    dim(spa) # 10557 structural variants x 37-1 genomes
    write_csv(spa, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/spa.csv"))

    # Make a longer list and Compute gene counts
    tt <- make_long_gpa(gpa); gpatl <- tt$gpatl; gene_order <- tt$gene_order
    dim(gpatl) # 244960 3
    dim(gene_order) # 26505 genes ordered by the prevalence across gneoms
    write_csv(gpatl, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/gpatl.csv"))
    write_csv(gene_order, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/gene_order.csv"))

    # Genes and the contigs they are from
    gd <- clean_gd(paste0(folder_data, "genomics/pangenome/", set_name,"/gene_data.csv"))
    dim(gd) # 258827 (gene x genome x contig) x 8 rows
    write_csv(gd, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/gd.csv"))

    # Full directory for each gene
    gpaf <- clean_gpaf(paste0(folder_data, "genomics/pangenome/", set_name,"/gene_presence_absence.csv"))
    dim(gpaf) # 26504 genes x 39-3 genomes. 1) gene name 2) non unique gene name 3) annotation
    gpafl <- make_long_gpaf(gpaf)
    dim(gpafl) # 244960 x 4

    # Make longer format of gene with contig info
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/gd.csv"))
    gpacl <- make_long_gpac(gpafl, contigs, gd)
    length(unique(gpacl$contig_id)) # 161 contigs
    write_csv(gpacl, paste0(folder_data, "genomics_analysis/gene_content/", set_name,"/gpacl.csv"))
}

# 1. All 36 genomes
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv")) %>% select(contig_id, replicon_type) %>% drop_na
clean_all("isolates", contigs)

# 2. 10 elevation medicae
clean_all("elev_med", contigs)

# 3. 17 urbanization meliloti
clean_all("urbn_mel", contigs)



