#' Compute gene-content divergence (GCV) as Dxy-like metric per replicon and region
#' Each gene = one site; alleles = presence (1) or absence (0)

library(tidyverse)
source(here::here("metadata.R"))

# paths ----
folder_ibd <- file.path(folder_genomics, "ibd/")

# input ----
gpacl <- read_csv(file.path(folder_ibd, "gpacl.csv"), show_col_types = FALSE)
isolates <- read_csv(file.path(folder_data, "mapping/isolates.csv"), show_col_types = FALSE)
ani <- read_csv(paste0(folder_genomics, "taxonomy/ani.csv"), show_col_types = FALSE)

# define groups ----
group_meliloti_PA <- isolates %>%
    left_join(ani) %>%
    filter(str_detect(organism_name, "meliloti"), str_detect(region, "Pennsylvania")) %>%
    filter(genome_id != "g42") %>%
    pull(genome_id)

group_medicae_VA <- isolates %>%
    left_join(ani) %>%
    filter(str_detect(organism_name, "medicae"), str_detect(region, "Virginia")) %>%
    pull(genome_id)

# prepare binary presence/absence matrix ----
gpa_replicon <- gpacl %>%
    filter(!is.na(replicon)) %>%
    mutate(present = 1) %>%
    distinct(gene, genome_id, replicon, present) %>%
    pivot_wider(names_from = genome_id, values_from = present, values_fill = 0)

# function to compute pairwise Dxy-like divergence for gene content ----
compute_gene_dxy <- function(df, genomes) {
    mat <- df %>% select(any_of(genomes)) %>% as.matrix()
    if (ncol(mat) < 2) return(NULL)

    # drop genes that are invariant (all 0 or all 1)
    variable_rows <- rowSums(mat) > 0 & rowSums(mat) < ncol(mat)
    mat <- mat[variable_rows, , drop = FALSE]
    if (nrow(mat) == 0) return(NULL)

    n_genes <- nrow(mat)
    combos <- combn(colnames(mat), 2, simplify = FALSE)

    tibble(
        genome_id1 = map_chr(combos, 1),
        genome_id2 = map_chr(combos, 2),
        gene_diff = map_dbl(combos, ~sum(mat[, .x[1]] != mat[, .x[2]])),
        total_genes = n_genes,
        Dxy_gene = gene_diff / n_genes
    )
}

# compute for each group and replicon ----
gcv_meliloti_PA <- gpa_replicon %>%
    group_split(replicon) %>%
    set_names(map_chr(., ~unique(.x$replicon))) %>%
    map_dfr(~compute_gene_dxy(.x, group_meliloti_PA), .id = "replicon") %>%
    mutate(group = "meliloti_PA")

gcv_medicae_VA <- gpa_replicon %>%
    group_split(replicon) %>%
    set_names(map_chr(., ~unique(.x$replicon))) %>%
    map_dfr(~compute_gene_dxy(.x, group_medicae_VA), .id = "replicon") %>%
    mutate(group = "medicae_VA")

# combine results ----
gcv <- bind_rows(gcv_meliloti_PA, gcv_medicae_VA)

# output ----
write_csv(gcv, file.path(folder_ibd, "gcv.csv"))
