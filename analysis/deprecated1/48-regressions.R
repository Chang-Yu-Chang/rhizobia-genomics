#' This script perform a couple regressions using the rhizobia site as resposne variable

library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
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
    mutate(clen_genome = clen_chromosome + clen_psyma + clen_psymb) %>%
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
    pivot_wider(id_cols = genome_name, names_from = contig_type, values_from = n_genes) %>%
    mutate(ngenes_genome = ngenes_chromosome + ngenes_psyma + ngenes_psymb)
# 0.3 clean the gene number for symbiosis genes
isolates_sym <- egcct %>%
    filter(str_detect(genome_name, "Chang")) %>%
    filter(str_detect(prokka_prodigal_acc, "nif|fix|nod")) %>%
    mutate(genome_name = factor(genome_name, paste0("Chang_Q5C_", 1:20))) %>%
    distinct(genome_name, contig_type, gene_cluster_id) %>%
    group_by(genome_name, contig_type) %>%
    count(name = "n_genes") %>%
    ungroup() %>%
    mutate(contig_type = paste0("nsymgenes_", contig_type)) %>%
    pivot_wider(id_cols = genome_name, names_from = contig_type, values_from = n_genes, values_fill = 0) %>%
    mutate(nsymgenes_genome = nsymgenes_chromosome + nsymgenes_psyma + nsymgenes_psymb)

# 0.3 join the four datasets
isolates <- isolates %>%
    left_join(isolates_egcc) %>%
    left_join(isolates_clen) %>%
    left_join(isolates_sym)


# 1. contig length and number of genes ----
p <- isolates %>%
    select(genome_name, strain_site_group, starts_with(c("ngenes", "clen"))) %>%
    pivot_longer(cols = -c(genome_name, strain_site_group), names_pattern = "(.*)_(.*)", names_to = c("name", "contig_type")) %>%
    pivot_wider(id_cols = c(genome_name, strain_site_group, contig_type)) %>%
    mutate(contig_type = factor(contig_type, c("genome", "chromosome", "psyma", "psymb"))) %>%
    ggplot() +
    geom_point(aes(x = clen/10^6, y = ngenes/1000, color = contig_type), stroke = 1, size = 3, shape = 21) +
    scale_color_aaas() +
    theme_classic() +
    theme(
        legend.position = c(0.2, 0.8),
        legend.background = element_rect(color = "black", fill = "white")
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = "contig length (Mbp)", y = "# of genes (x1000)")

ggsave(paste0(folder_data, "temp/48-01-contig_size.png"), p, width = 4, height = 4)

# 2. number of total genes vs growth trait ----
# Does the total number of genes in a genome predict its growth traits?
isolates_ngenes <- isolates %>%
    select(genome_id, r, lag, maxOD, starts_with(c("ngenes"))) %>%
    pivot_longer(cols = starts_with("ngenes"), names_prefix = "ngenes_", names_to = "contig_type", values_to = "ngenes") %>%
    pivot_longer(cols = c(r, lag, maxOD), names_to = "growth_trait", values_to = "value")
p <- isolates_ngenes %>%
    ggplot() +
    geom_point(aes(x = ngenes, y = value), stroke = 1, size = 3, shape = 21) +
    geom_smooth(aes(x = ngenes, y = value), method = "lm") +
    facet_grid(growth_trait ~ contig_type, scales = "free") +
    #scale_color_aaas() +
    theme_bw() +
    theme(
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank()
    ) +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/48-02-ngenes_vs_r.png"), p, width = 12, height = 9)

# Use purrr to fit regression models for each group
fit_regression <- function(data) {
    model <- lm(value ~ ngenes, data = data)
    return(broom::tidy(model))
}

isolates_ngenes %>%
    filter(contig_type == "chromosome", growth_trait == "r") %>%
    select(genome_id, ngenes, value) %>%
    fit_regression()

isolates_ngenes %>%
    nest(data = c(genome_id, ngenes, value)) %>%
    mutate(coefficient = map(data, fit_regression)) %>%
    select(contig_type, growth_trait, coefficient) %>%
    unnest(c(coefficient)) %>%
    filter(term == "ngenes") %>%
    filter(p.value < 0.05) %>%
    select(contig_type, growth_trait, estimate, p.value)

# 3. Does the contig size predict its growth? ----
fit_regression <- function(data) {
    model <- lm(value ~ clen, data = data)
    return(broom::tidy(model))
}


isolates_clen <- isolates %>%
    select(genome_id, r, lag, maxOD, starts_with(c("clen_"))) %>%
    pivot_longer(cols = starts_with("clen_"), names_prefix = "clen_", names_to = "contig_type", values_to = "clen") %>%
    pivot_longer(cols = c(r, lag, maxOD), names_to = "growth_trait", values_to = "value")

isolates_clen %>%
    filter(contig_type == "chromosome", growth_trait == "r") %>%
    select(genome_id, clen, value) %>%
    fit_regression()

isolates_clen %>%
    nest(data = c(genome_id, clen, value)) %>%
    mutate(coefficient = map(data, fit_regression)) %>%
    select(contig_type, growth_trait, coefficient) %>%
    unnest(c(coefficient)) %>%
    filter(term == "clen") %>%
    filter(p.value < 0.05) %>%
    select(contig_type, growth_trait, estimate, p.value)


# 4. Does the number of symbiosis gene predict its growth? ----
isolates_nsymgenes <- isolates %>%
    select(genome_id, r, lag, maxOD, starts_with(c("nsymgenes"))) %>%
    pivot_longer(cols = starts_with("nsymgenes"), names_prefix = "nsymgenes_", names_to = "contig_type", values_to = "nsymgenes") %>%
    pivot_longer(cols = c(r, lag, maxOD), names_to = "growth_trait", values_to = "value")
p <- isolates_nsymgenes %>%
    ggplot() +
    geom_point(aes(x = nsymgenes, y = value), stroke = 1, size = 3, shape = 21) +
    geom_smooth(aes(x = nsymgenes, y = value), method = "lm") +
    facet_grid(growth_trait ~ contig_type, scales = "free") +
    #scale_color_aaas() +
    theme_bw() +
    theme(
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank()
    ) +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/48-04-nsymgenes_vs_r.png"), p, width = 12, height = 9)

# Fit logistic regression model
fit_regression <- function(data) {
    model <- lm(value ~ nsymgenes, data = data)
    return(broom::tidy(model))
}

isolates_nsymgenes %>%
    filter(contig_type == "chromosome", growth_trait == "r") %>%
    select(genome_id, nsymgenes, value) %>%
    fit_regression()

isolates_nsymgenes %>%
    nest(data = c(genome_id, nsymgenes, value)) %>%
    mutate(coefficient = map(data, fit_regression)) %>%
    select(contig_type, growth_trait, coefficient) %>%
    unnest(c(coefficient)) %>%
    filter(term == "nsymgenes") %>%
    filter(p.value < 0.05) %>%
    select(contig_type, growth_trait, estimate, p.value)



























