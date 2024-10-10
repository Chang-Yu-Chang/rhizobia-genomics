#' This script plots the heatmap of gene content

library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))

    return(list(gpa = gpa, gene_order = gene_order, gpatl = gpatl, sml = sml))
}
plot_heatmap <- function (gpatl) {
    gpatl %>%
        ggplot() +
        geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
        scale_y_discrete(expand = c(0,0)) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 5, color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        ) +
        guides(fill = "none") +
        labs(x = "gene cluster", y = "genome")
}
plot_gfs <- function (gpa) {
    #' Plot the gene frequency spectrum
    gfs <- gpa %>%
        pivot_longer(-gene) %>%
        group_by(gene) %>%
        summarize(n_genomes = sum(value)) %>%
        group_by(n_genomes) %>%
        count()
    gfs %>%
        ggplot() +
        geom_col(aes(x = n_genomes, y = n), color = "black", fill = "white") +
        scale_x_continuous(breaks = c(1, 31, seq(5, 30, 5))) +
        theme_classic() +
        theme() +
        guides() +
        labs(x = "# of genomes", y = "# of genes")
}
plot_singleton <- function (gpa, gpatl, isolates_tax) {
    #' Plot the number of singletons per genome
    list_sg <- gpa$gene[which(apply(gpa[,-1], 1, sum) == 1)]
    gpatlsg <- gpatl %>%
        filter(gene %in% list_sg) %>%
        left_join(isolates_tax) %>%
        mutate(contig_species = str_remove(contig_species, "E. ")) %>%
        mutate(genome_id = factor(genome_id, isolates$genome_id))


    gpatlsg %>%
        group_by(genome_id, contig_species) %>%
        count() %>%
        ggplot() +
        geom_col(aes(x = genome_id, y = n, fill = contig_species), color = "black") +
        scale_fill_manual(values = species_colors) +
        theme_bw() +
        theme() +
        guides() +
        labs(title = "Symbiotic and non-symbiotic strains", y = "# of singltons")
}
plot_bcsm <- function (sml) {
    #' Plot the bray-curtis similarity matrix. Input is in long format
    sml %>%
        mutate(genome_id1 = factor(genome_id1, isolates_tax$genome_id), genome_id2 = factor(genome_id2, rev(isolates_tax$genome_id))) %>%
        ggplot() +
        geom_tile(aes(x = genome_id1, y = genome_id2, fill = bray_curtis_similarity)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradient(low = "snow", high = "maroon") +
        theme_bw() +
        theme() +
        labs()
}

# 1. 36 isolate genomes ----
tt <- read_gpas("isolates")
nrow(tt$gene_order)
p <- plot_heatmap(tt$gpatl)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/isolates-01-gpa_heatmap.png"), p, width = 6, height = 3)
p <- plot_gfs(tt$gpa)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/isolates-02-gene_frequency_spectrum.png"), p, width = 5, height = 4)
p <- plot_singleton(tt$gpa, tt$gpatl, isolates_tax)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/isolates-03-singletons.png"), p, width = 10, height = 8)
p <- plot_bcsm(tt$sml)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/isolates-04-gpa_sm.png"), p, width = 10, height = 8)


# 2. Elevation medicae  ----
tt <- read_gpas("elev_med")
nrow(tt$gene_order)
p <- plot_heatmap(tt$gpatl)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/elev_med-01-gpa_heatmap.png"), p, width = 6, height = 3)
p <- plot_gfs(tt$gpa)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/elev_med-02-gene_frequency_spectrum.png"), p, width = 5, height = 4)
p <- plot_singleton(tt$gpa, tt$gpatl, isolates_tax)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/elev_med-03-singletons.png"), p, width = 10, height = 8)
p <- plot_bcsm(tt$sml)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/elev_med-04-gpa_sm.png"), p, width = 10, height = 8)

# 3. Urbanization meliloti  ----
tt <- read_gpas("urbn_mel")
nrow(tt$gene_order)
p <- plot_heatmap(tt$gpatl)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/urbn_mel-01-gpa_heatmap.png"), p, width = 6, height = 3)
p <- plot_gfs(tt$gpa)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/urbn_mel-02-gene_frequency_spectrum.png"), p, width = 5, height = 4)
p <- plot_singleton(tt$gpa, tt$gpatl, isolates_tax)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/urbn_mel-03-singletons.png"), p, width = 10, height = 8)
p <- plot_bcsm(tt$sml)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/urbn_mel-04-gpa_sm.png"), p, width = 10, height = 8)
























