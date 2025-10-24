#' Compute per-replicon Dxy using stable SCCG and population-level SNP sites

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

# paths ----
folder_core <- file.path(folder_genomics, "pangenome/panaroo/aligned_gene_sequences/")
folder_ibd  <- file.path(folder_genomics, "ibd/")

# input ----
list_sccg <- read_csv(paste0(folder_genomics, "pangenome/gene_content/list_sccg.csv"), col_names = "gene", show_col_types = FALSE)
gpacl     <- read_csv(file.path(folder_ibd, "gpacl.csv"), show_col_types = FALSE)
isolates  <- read_csv(file.path(folder_data, "mapping/isolates.csv"), show_col_types = FALSE)
ani       <- read_csv(paste0(folder_genomics, "taxonomy/ani.csv"), show_col_types = FALSE)

# define groups ----
group_meliloti_PA <- isolates %>%
    left_join(ani) %>%
    filter(str_detect(organism_name, "meliloti"), str_detect(region, "Pennsylvania"), genome_id != "g42") %>%
    pull(genome_id)

group_medicae_VA <- isolates %>%
    left_join(ani) %>%
    filter(str_detect(organism_name, "medicae"), str_detect(region, "Virginia")) %>%
    pull(genome_id)

# identify stable-replicon SCCGs ----
get_stable_sccg <- function(genomes) {
    replicon_map <- gpacl %>%
        filter(gene %in% list_sccg$gene, genome_id %in% genomes, !is.na(replicon)) %>%
        distinct(gene, genome_id, replicon)

    stable_sccg <- replicon_map %>%
        group_by(gene) %>%
        summarise(n_replicon = n_distinct(replicon), .groups = "drop") %>%
        filter(n_replicon == 1) %>%
        pull(gene)

    replicon_map %>%
        filter(gene %in% stable_sccg) %>%
        distinct(gene, replicon)
}

stable_meliloti_PA <- get_stable_sccg(group_meliloti_PA)
stable_medicae_VA  <- get_stable_sccg(group_medicae_VA)

# helper: count SNP sites across population ----
count_population_snps <- function(file, genomes) {
    aln <- tryCatch(read.dna(file, format = "fasta"), error = function(e) NULL)
    if (is.null(aln)) return(NULL)

    new_names <- rownames(aln) %>% str_extract("^[^;]+")
    rownames(aln) <- new_names
    genomes_present <- intersect(rownames(aln), genomes)
    if (length(genomes_present) < 2) return(NULL)

    aln <- aln[genomes_present, , drop = FALSE]
    if (ncol(aln) == 0) return(NULL)

    m <- as.matrix(aln)
    variable_sites <- apply(m, 2, function(col) {
        nuc <- col[col != "-"]
        length(unique(nuc)) > 1
    })

    tibble(
        gene = tools::file_path_sans_ext(basename(file)),
        n_sites_total = ncol(m),
        n_sites_poly  = sum(variable_sites)
    ) %>%
        mutate(gene = str_remove(gene, ".aln"))
}

# helper: compute pairwise SNP Dxy ----
compute_snp_dxy <- function(file, genomes) {
    aln <- tryCatch(read.dna(file, format = "fasta"), error = function(e) NULL)
    if (is.null(aln)) return(NULL)

    new_names <- rownames(aln) %>% str_extract("^[^;]+")
    rownames(aln) <- new_names
    genomes_present <- intersect(rownames(aln), genomes)
    if (length(genomes_present) < 2) return(NULL)

    aln <- aln[genomes_present, , drop = FALSE]
    m <- as.matrix(aln)
    d <- dist.dna(m, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)

    tibble(
        gene = tools::file_path_sans_ext(basename(file)),
        genome_id1 = rep(rownames(d), each = nrow(d)),
        genome_id2 = rep(rownames(d), times = nrow(d)),
        snp_diff = as.vector(d)
    ) %>%
        mutate(gene = str_remove(gene, ".aln")) %>%
        filter(genome_id1 < genome_id2)
}

# list alignments ----
alignments <- list.files(folder_core, pattern = "\\.fas$", full.names = TRUE)

# function to run Dxy per group ----
process_group <- function(genomes, stable_list, tag) {
    files <- alignments[basename(alignments) %in% paste0(stable_list$gene, ".aln.fas")]

    snp_table <- map_dfr(files, compute_snp_dxy, genomes = genomes)
    snp_sites <- map_dfr(files, count_population_snps, genomes = genomes)

    # total number of polymorphoric sites
    total_sites <- snp_sites %>%
        left_join(stable_list) %>%
        group_by(replicon) %>%
        summarize(total_sites_poly = sum(n_sites_poly))

    snp_table %>%
        left_join(stable_list, by = "gene") %>%
        group_by(replicon, genome_id1, genome_id2) %>%
        summarise(total_snp_diff = sum(snp_diff, na.rm = TRUE), .groups = "drop") %>%
        left_join(total_sites) %>%
        mutate(
            Dxy_snp = total_snp_diff / total_sites_poly,
            group = tag
        )
}

# run for both groups ----
sccg_meliloti_PA <- process_group(group_meliloti_PA, stable_meliloti_PA, "meliloti_PA")
sccg_medicae_VA  <- process_group(group_medicae_VA,  stable_medicae_VA,  "medicae_VA")

# combine ----
sccg <- bind_rows(sccg_meliloti_PA, sccg_medicae_VA)

# output ----
write_csv(sccg, file.path(folder_ibd, "sccg.csv"))
