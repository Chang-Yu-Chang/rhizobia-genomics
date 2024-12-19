#' Compute trees: gpa,

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
compute_gpa_tree <- function (tb) {
    d <- dist(t(tb[,-1]))
    tr <- hclust(d, method = "ward.D2") %>% as.phylo()
    return(tr)
}

n_genes <- list()
for (set_name in c("elev_med", "urbn_mel")) {
    tt <- read_gpas(set_name)
    gpac <- tt$gpacl %>%
        select(replicon_type, genome_id, gene) %>%
        group_by(replicon_type) %>%
        mutate(temp = 1) %>%
        pivot_wider(names_from = genome_id, values_from = temp, values_fill = 0) %>%
        drop_na(replicon_type) %>%
        nest(data = -replicon_type) %>%
        mutate(tr = map(data, compute_gpa_tree))

    for (replicon in c("chromosome", "pSymA", "pSymB", "pAcce")) {
        tr_gpa <- gpac$tr[[which(gpac$replicon_type == replicon)]]
        write.tree(tr_gpa, paste0(folder_data, "phylogenomics_analysis/replicon_trees/", set_name, "/", replicon, "/tr_gpa.tree"))
    }

    count_n_core <- function (y) sum((apply(y[,-1], 1, function (x) all(x==1))))
    count_n_accessory <- function (y) sum((apply(y[,-1], 1, function (x) !all(x==1))))
    n_genes[[set_name]] <- gpac %>% mutate(
        n_core = map_int(data, count_n_core),
        n_accessory = map_int(data, count_n_accessory),
    )
}

# Number of core vs accessory genes on
n_genes <- n_genes %>%
    bind_rows(.id = "set_name") %>%
    select(set_name, replicon_type, n_core, n_accessory)
write_csv(n_genes, paste0(folder_data, "phylogenomics_analysis/replicon_trees/n_accessory.csv"))
