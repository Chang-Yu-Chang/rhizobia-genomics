#' Compute Mantel tests for core (SCCG) and accessory (GCV) Dxy and IBD

library(tidyverse)
library(vegan)
source(here::here("metadata.R"))

# input ----
sccg <- read_csv(file.path(folder_genomics, "ibd/sccg.csv"))
gcv  <- read_csv(file.path(folder_genomics, "ibd/gcv.csv"))
isolates <- read_csv(file.path(folder_data, "mapping/isolates.csv"))
site_dist <- read_csv(file.path(folder_phenotypes, "sites/sites_dist.csv"))  # columns: genome_id1, genome_id2, dist_km

# combine datasets ----
sccg_long <- sccg %>%
    rename(Dxy = Dxy_snp) %>%
    mutate(genome_pair = paste(replicon, pmin(genome_id1, genome_id2), pmax(genome_id1, genome_id2), sep = "_"))

gcv_long <- gcv %>%
    rename(Dxy = Dxy_gene) %>%
    mutate(genome_pair = paste(replicon, pmin(genome_id1, genome_id2), pmax(genome_id1, genome_id2), sep = "_"))

dist_full <- crossing(genome_id1 = isolates$genome_id, genome_id2 = isolates$genome_id) %>%
    left_join(isolates %>% select(genome_id1 = genome_id, site1 = site)) %>%
    left_join(isolates %>% select(genome_id2 = genome_id, site2 = site)) %>%
    left_join(site_dist)

dist_long <- select(gcv_long, replicon, genome_id1, genome_id2) %>%
    left_join(dist_full) %>%
    mutate(genome_pair = paste(replicon, pmin(genome_id1, genome_id2), pmax(genome_id1, genome_id2), sep = "_")) %>%
    select(replicon, genome_pair, dist_geo_km)


# join ----
sccg_join <- sccg_long %>% left_join(dist_long)
gcv_join <- gcv_long %>% left_join(dist_long)

# function for Mantel test ----
run_mantel <- function(df) {
    if (nrow(df) < 3) return(NULL)

    df <- df %>%
        distinct(genome_id1, genome_id2, Dxy, dist_geo_km, .keep_all = TRUE)

    genomes <- unique(c(df$genome_id1, df$genome_id2))
    mat_dxy <- matrix(NA, nrow = length(genomes), ncol = length(genomes),
                      dimnames = list(genomes, genomes))
    mat_dist <- mat_dxy

    for (i in seq_len(nrow(df))) {
        g1 <- df$genome_id1[i]; g2 <- df$genome_id2[i]
        mat_dxy[g1, g2] <- mat_dxy[g2, g1] <- df$Dxy[i]
        mat_dist[g1, g2] <- mat_dist[g2, g1] <- df$dist_geo_km[i]
    }

    res <- vegan::mantel(mat_dxy, mat_dist, method = "pearson", permutations = 999)
    tibble(r = res$statistic, p = res$signif)
}

# summarize Mantel per group × replicon ----
mantel_sccg <- sccg_join %>%
    group_by(group, replicon) %>%
    group_modify(~run_mantel(.x)) %>%
    mutate(type = "Core genome (sequence similarity)")

mantel_gcv <- gcv_join %>%
    group_by(group, replicon) %>%
    group_modify(~run_mantel(.x)) %>%
    mutate(type = "Accessory genome (gene content similarity)")

mantel_results <- bind_rows(mantel_sccg, mantel_gcv)
write_csv(sccg_join, file.path(folder_genomics, "ibd/sccg_join.csv"))
write_csv(gcv_join, file.path(folder_genomics, "ibd/gcv_join.csv"))
write_csv(mantel_results, file.path(folder_genomics, "ibd/mantel_results.csv"))
