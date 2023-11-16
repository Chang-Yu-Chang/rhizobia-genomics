#' This script perform a couple regressions using the rhizobia site as resposne variable

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))


# 0. read data ----

# contig length
g_contigs <- read_csv(paste0(folder_data, "temp/34-g_contigs.csv"), show_col_types = F)
# contig data and gene numbers
egcct <- read_csv(paste0(folder_data, "temp/45-egcct.csv"), show_col_types = F)
# isolate traits
isolates <- read_csv(paste0(folder_data, "temp/42-isolates.csv"), show_col_types = F) %>%
    mutate(genome_name = str_replace(genome_id, "g", "Chang_Q5C_") %>% factor(paste0("Chang_Q5C_", 1:20)))

# 0.1 clean the contig length
t_contigs <- egcct %>% distinct(genome_name, contig, contig_type) # contigtypes for cross reference
isolates_clen <- g_contigs %>%
    mutate(genome_name = str_replace(genome_id, "g", "Chang_Q5C_") %>% factor(paste0("Chang_Q5C_", 1:20))) %>%
    left_join(t_contigs) %>%
    drop_na %>% # drop those contigs that do not match to the three elements. Drop the ncbi genomes too
    mutate(genome_name = factor(genome_name, paste0("Chang_Q5C_", 1:20))) %>%
    select(genome_name, contig_type, contig_length) %>%
    mutate(contig_type = paste0("clen_", contig_type)) %>%
    pivot_wider(id_cols = genome_name, names_from = contig_type, values_from = contig_length) %>%
    arrange(genome_name)
# 0.2 Clean the gene number
isolates_egcc <- egcct %>%
    filter(str_detect(genome_name, "Chang")) %>%
    mutate(genome_name = factor(genome_name, paste0("Chang_Q5C_", 1:20))) %>%
    distinct(genome_name, contig_type, gene_cluster_id) %>%
    group_by(genome_name, contig_type) %>%
    count(name = "n_genes") %>%
    ungroup() %>%
    mutate(contig_type = paste0("ngenes_", contig_type)) %>%
    pivot_wider(id_cols = genome_name, names_from = contig_type, values_from = n_genes)
# 0.3 join the three datasets
isolates <- isolates %>%
    left_join(isolates_egcc) %>%
    left_join(isolates_clen)


# 1.
p <- isolates %>%
    select(genome_name, strain_site_group, starts_with(c("ngenes", "clen"))) %>%
    pivot_longer(cols = -c(genome_name, strain_site_group), names_pattern = "(.*)_(.*)", names_to = c("name", "contig_type")) %>%
    pivot_wider(id_cols = c(genome_name, strain_site_group, contig_type)) %>%
    ggplot() +
    geom_point(aes(x = clen/10^6, y = ngenes/1000, color = contig_type), stroke = 1, size = 3, shape = 21) +
    scale_color_aaas() +
    theme_classic() +
    theme(
        legend.position = c(0.2, 0.8),
        legend.background = element_rect(color = "black", fill = "white")
    ) +
    guides() +
    labs(x = "contig length (Mbp)", y = "# of genes (x1000)")

ggsave(paste0(folder_data, "temp/48-01-contig_size.png"), p, width = 4, height = 4)











