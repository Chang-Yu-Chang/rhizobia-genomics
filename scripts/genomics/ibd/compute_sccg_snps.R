#' Compute SNP-based Dxy per replicon and region
#' using only stable-replicon single-copy core genes (SCCGs)
#' Mirrors the GCV script

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

# paths ----
folder_core <- file.path(folder_genomics, "pangenome/panaroo/aligned_gene_sequences/")
folder_ibd <- file.path(folder_genomics, "ibd/")

# input ----
list_sccg <- read_csv(paste0(folder_genomics, "pangenome/gene_content/list_sccg.csv"),
                      col_names = "gene", show_col_types = FALSE)
gpacl <- read_csv(file.path(folder_ibd, "gpacl.csv"), show_col_types = FALSE)
isolates <- read_csv(file.path(folder_data, "mapping/isolates.csv"), show_col_types = FALSE)
ani <- read_csv(paste0(folder_genomics, "taxonomy/ani.csv"), show_col_types = FALSE)

# define groups ----
group_meliloti_PA <- isolates %>%
    left_join(ani) %>%
    filter(str_detect(organism_name, "meliloti"), str_detect(region, "Pennsylvania")) %>%
    pull(genome_id)

group_medicae_VA <- isolates %>%
    left_join(ani) %>%
    filter(str_detect(organism_name, "medicae"), str_detect(region, "Virginia")) %>%
    pull(genome_id)

# identify stable-replicon SCCGs ----
replicon_map <- gpacl %>%
    filter(gene %in% list_sccg$gene, !is.na(replicon)) %>%
    distinct(gene, genome_id, replicon)

stable_sccg <- replicon_map %>%
    group_by(gene) %>%
    summarise(n_replicon = n_distinct(replicon), .groups = "drop") %>%
    filter(n_replicon == 1) %>%
    pull(gene)

list_sccg_stable <- replicon_map %>%
    filter(gene %in% stable_sccg) %>%
    distinct(gene, replicon)

write_csv(list_sccg_stable, file.path(folder_ibd, "list_sccg_stable.csv"))

# function to compute per-gene pairwise SNP distances ----
compute_snp_dxy <- function(file, genomes) {
    aln <- tryCatch(read.dna(file, format = "fasta"), error = function(e) NULL)
    if (is.null(aln)) return(NULL)

    new_names <- rownames(aln) %>% str_extract("^[^;]+")
    rownames(aln) <- new_names

    genomes_present <- intersect(rownames(aln), genomes)
    if (length(genomes_present) < 2) return(NULL)

    aln <- aln[genomes_present, , drop = FALSE]
    if (ncol(aln) == 0) return(NULL)

    m <- as.matrix(aln)
    d <- dist.dna(m, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)

    tibble(
        gene = tools::file_path_sans_ext(basename(file)),
        genome_id1 = rep(rownames(d), each = nrow(d)),
        genome_id2 = rep(rownames(d), times = nrow(d)),
        snp_diff = as.vector(d),
        n_sites = ncol(m),
        Dxy_snp = snp_diff / n_sites
    ) %>%
        filter(genome_id1 < genome_id2)
}

# list SCCG alignments ----
alignments <- list.files(folder_core, pattern = "\\.fas$", full.names = TRUE)
sccg_alignments <- alignments[basename(alignments) %in% paste0(stable_sccg, ".aln.fas")]

# compute pairwise SNP Dxy ----
sccg_meliloti_PA <- map_dfr(sccg_alignments, compute_snp_dxy, genomes = group_meliloti_PA) %>%
    mutate(group = "meliloti_PA")

sccg_medicae_VA <- map_dfr(sccg_alignments, compute_snp_dxy, genomes = group_medicae_VA) %>%
    mutate(group = "medicae_VA")

# combine and join replicon info ----
sccg <- bind_rows(sccg_meliloti_PA, sccg_medicae_VA) %>%
    mutate(gene = str_remove(gene, ".aln")) %>%
    left_join(list_sccg_stable, by = "gene")

# aggregate per-pair per-replicon ----
sccg_ibd <- sccg %>%
    group_by(group, replicon, genome_id1, genome_id2) %>%
    summarise(
        total_snp_diff = sum(snp_diff, na.rm = TRUE),
        total_sites = sum(n_sites, na.rm = TRUE),
        Dxy_snp = total_snp_diff / total_sites,
        .groups = "drop"
    )

# output ----
write_csv(sccg, file.path(folder_ibd, "sccg.csv"))
write_csv(sccg_ibd, file.path(folder_ibd, "sccg_dxy.csv"))
