#'

library(tidyverse)
library(janitor)
library(cowplot)
library(vegan)
library(ggh4x)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
symbiosis_genes <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/symbiosis_genes.csv"))
tt <- read_gpas()


tb <- tibble(genome_id = factor(iso$genome_id[1:23], iso$genome_id)) %>%
    left_join(iso) %>%
    mutate(cog_classify = map(genome_id, ~read_tsv(paste0(folder_data, "genomics/cog/", .x, "/cog_classify.tsv"))))


tb %>%
    unnest(cog_classify) %>%
    clean_names() %>%
    group_by(contig_species, genome_id, cog_letter) %>%
    #filter(cog_id == "COG0071") %>%
    filter(str_detect(cog_name, " heat ")) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = n, fill = cog_letter), color = 1) +
    facet_grid(~contig_species, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()


