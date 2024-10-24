#' This script plots the heatmap of gene content

library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

make_acce_fst <- function (per_acce_fst, gene_order, gpacl) {

    gene_replicon <- filter(gpacl, genome_id == gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))
    per_acce_fst %>%
        mutate(gene = factor(gene, tt$gene_order$gene)) %>%
        left_join(gene_replicon) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce")))
}
plot_heatmap <- function (gpa, gpatl, gpacl, list_wgpa, by_replicon = F) {
    list_cg <- gpa$gene[apply(gpa[,-1], 1, sum) == ncol(gpa)-1]
    gene_order <- unique(gpacl$gene)
    p <- gpacl %>%
        filter(!gene %in% list_cg) %>%
        # Exclude those genes that are assigned to different contigs in differnt genomes
        #filter(!gene %in% list_wgpa$gene) %>%
        left_join(isolates) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
        ggplot() +
        geom_tile(aes(x = gene, y = genome_id, fill = population)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_manual(values = population_colors) +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 10, color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
        ) +
        guides() +
        labs(x = "accessory gene", y = "genome")

    if (by_replicon) {
        return(p + facet_grid(.~replicon_type, scales = "free_x", space = "free_x"))
    } else return(p)

}
plot_acce_fst <- function (acce_fst) {
    #' Plot accessory gene fst
    acce_fst %>%
        ggplot() +
        geom_segment(aes(x = gene, xend = gene, y = fst-0.01, yend = fst+0.01), alpha = 0.4) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
        ) +
        guides() +
        labs(x = "accessory gene", y = "Fst")
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
plot_singleton <- function (gpa, gpatl) {
    #' Plot the number of singletons per genome
    list_sg <- gpa$gene[which(apply(gpa[,-1], 1, sum) == 1)]
    gpatlsg <- gpatl %>%
        filter(gene %in% list_sg) %>%
        mutate(genome_id = factor(genome_id, isolates$genome_id))

    gpatlsg %>%
        group_by(genome_id) %>%
        count() %>%
        ggplot() +
        geom_col(aes(x = genome_id, y = n), color = "black") +
        theme_bw() +
        theme() +
        guides() +
        labs(title = "Symbiotic and non-symbiotic strains", y = "# of singltons")
}
plot_bcsm <- function (sml) {
    #' Plot the bray-curtis similarity matrix. Input is in long format
    sml %>%
        mutate(genome_id1 = factor(genome_id1, isolates$genome_id), genome_id2 = factor(genome_id2, rev(isolates$genome_id))) %>%
        ggplot() +
        geom_tile(aes(x = genome_id1, y = genome_id2, fill = bray_curtis_similarity)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradient(low = "snow", high = "maroon") +
        theme_bw() +
        theme() +
        labs()
}

for (set_name in c("elev_med", "urbn_mel")) {
    #set_name = "elev_med"
    #set_name = "urbn_mel"
    tt <- read_gpas(set_name)
    per_acce_fst <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/per_acce_fst.csv"))
    n_all <- nrow(tt$gene_order) # number of all genes
    n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes
    acce_fst <- make_acce_fst(per_acce_fst, tt$gene_order, tt$gpacl)


    p1 <- plot_heatmap(tt$gpa, tt$gpatl, tt$gpacl, list_wgpa) + ggtitle(paste0(n_accessory, " accessory genes"))
    p2 <- plot_acce_fst(acce_fst)
    p <- plot_grid(p1, p2, nrow = 2, axis = "lrt", align = "vh")
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-01-gpa_heatmap.png"), p, width = 10, height = 6)
    p <- plot_gfs(tt$gpa)
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-02-gene_frequency_spectrum.png"), p, width = 5, height = 4)
    p <- plot_singleton(tt$gpa, tt$gpatl)
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-03-singletons.png"), p, width = 10, height = 8)
    p <- plot_bcsm(tt$sml)
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-04-gpa_sm.png"), p, width = 10, height = 8)
}

if (F) {

# Genes that are assigned to different contigs in different genomes
list_wgpa <- tt$gpacl %>%
    select(gene, genome_id, replicon_type) %>%
    distinct(gene, genome_id, .keep_all = T) %>%
    group_by(gene) %>%
    summarize(is_all_same = length(unique(replicon_type[!is.na(replicon_type)])) %in% c(0,1)) %>%
    filter(!is_all_same)

}

