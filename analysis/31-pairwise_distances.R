#' This script joins the distances

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
dist_genetics <- read_csv(paste0(folder_data, 'temp/19-dist_genetics.csv'))
dist_traits <- read_csv(paste0(folder_data, "temp/29-dist_traits.csv"))

# Join the distances in genetics and traits
dists <- dist_genetics %>% left_join(dist_traits)
dists_long <- dists %>%
    pivot_longer(cols = c(-genome_id1, -genome_id2), names_to = "d_type", names_prefix = "d_") %>%
    mutate(d_group = case_when(
        d_type %in% c("ani", "kmer", "jaccard", "fluidity") ~ "genetic",
        d_type %in% c("growth", "symbiosis") ~ "composite trait",
        d_type == "geo" ~ "geographic",
        TRUE ~ "trait" 
    )) %>%
    select(genome_id1, genome_id2, d_group, d_type, value)


write_csv(dists, paste0(folder_data, "temp/31-dists.csv"))
write_csv(dists_long, paste0(folder_data, "temp/31-dists_long.csv"))

# Correlation between traits
list_genomics <- c("ani", "kmer", "jaccard", "fluidity")
list_traits <- c("dry_weight_mg", "nodule_number", "root_weight_mg",
        str_subset(colnames(isolates_traits), "r_"),
        str_subset(colnames(isolates_traits), "lag_"),
        str_subset(colnames(isolates_traits), "maxOD_"))

tb_gt <- crossing(d_g = list_genomics, d_t = list_traits) %>%
    rowwise %>%
    mutate(dat = list(tibble(dists[,paste0("d_", c(d_g, d_t))]))) %>%
    ungroup()


compute_cor <- function(tb) {
    #' The input is a two-column tibble
    t1 <- unlist(tb[,1])
    t2 <- unlist(tb[,2])
    result <- cor.test(t1, t2, method = "pearson") %>% broom::tidy()
    return(result)
}

tb_gt <- tb_gt %>%
    rowwise() %>%
    mutate(result = list(compute_cor(dat))) %>%
    select(-dat) %>%
    unnest(cols = result)
write_csv(tb_gt, paste0(folder_data, "temp/31-tb_gt.csv"))
