#' Assign the gpa file to replicon based on which contig its on

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

# paths ----
folder_core <- file.path(folder_genomics, "pangenome/panaroo/aligned_gene_sequences")
folder_ibd <- file.path(folder_genomics, "ibd")

# input ----
list_sccg <- read_csv(paste0(folder_genomics, "pangenome/gene_content/list_sccg.csv"), col_names = "gene", show_col_types = FALSE)
gpacl <- read_csv(file.path(folder_ibd, "gpacl.csv"), show_col_types = FALSE)
isolates <- read_csv(file.path(folder_genomics, "mapping/isolates.csv"), show_col_types = FALSE)

# define groups ----
group_meliloti_PA <- isolates %>%
    filter(str_detect(organism_name, "meliloti"), str_detect(region, "Pennsylvania")) %>%
    pull(genome_id)

group_medicae_VA <- isolates %>%
    filter(str_detect(organism_name, "medicae"), str_detect(region, "Virginia")) %>%
    pull(genome_id)


# function to count snps ----
count_snps <- function(file, genomes) {
    aln <- tryCatch(read.dna(file, format = "fasta"), error = function(e) NULL)
    if (is.null(aln)) return(NULL)
    genomes_present <- intersect(rownames(aln), genomes)
    if (length(genomes_present) < 2) return(NULL)
    aln <- aln[genomes_present, , drop = FALSE]
    if (ncol(aln) == 0) return(NULL)

    m <- as.matrix(aln)
    n_sites <- ncol(m)
    variable_sites <- apply(m, 2, function(col) {
        nuc <- col[col != "-"]
        length(unique(nuc)) > 1
    })

    tibble(
        gene = tools::file_path_sans_ext(basename(file)),
        n_sites = n_sites,
        snp_count = sum(variable_sites)
    )
}

# list core alignments ----
alignments <- list.files(folder_core, pattern = "\\.fa(sta)?$", full.names = TRUE)
sccg_alignments <- alignments[basename(alignments) %in% paste0(list_sccg$gene, ".fasta")]

# compute snps for each region ----
snp_meliloti_PA <- map_dfr(sccg_alignments, count_snps, genomes = group_meliloti_PA) %>%
    mutate(group = "meliloti_PA")

snp_medicae_VA <- map_dfr(sccg_alignments, count_snps, genomes = group_medicae_VA) %>%
    mutate(group = "medicae_VA")

snp_table <- bind_rows(snp_meliloti_PA, snp_medicae_VA)

# assign replicon ----
snp_table <- snp_table %>%
    left_join(gpacl %>% select(gene, replicon) %>% distinct(), by = "gene") %>%
    filter(!is.na(replicon))

# summarize by group and replicon ----
snp_summary <- snp_table %>%
    group_by(group, replicon) %>%
    summarise(
        n_genes = n(),
        total_snps = sum(snp_count, na.rm = TRUE),
        total_sites = sum(n_sites, na.rm = TRUE),
        mean_snps_per_gene = mean(snp_count, na.rm = TRUE),
        snp_rate = total_snps / total_sites,
        .groups = "drop"
    )
