#' This script plots the heatmap of gene content

library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- arrange(isolates, site_group)

read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/spa.csv"))

    return(list(gpa = gpa, gene_order = gene_order, gpatl = gpatl, gpacl = gpacl, gd = gd, sml = sml, list_sccg = list_sccg, spa = spa))
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
        geom_tile(aes(x = gene, y = genome_id, fill = site_group)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_manual(values = site_group_colors) +
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
plot_fisher <- function (tidied_fisher, gpa, gpatl, gpacl, list_wgpa, by_replicon = F) {
    list_cg <- gpa$gene[apply(gpa[,-1], 1, sum) == ncol(gpa)-1]
    gene_order <- unique(gpatl$gene)
    p <- tidied_fisher %>%
        filter(!gene %in% list_cg) %>%
        # Exclude those genes that are assigned to different contigs in differnt genomes
        #filter(!gene %in% list_wgpa$gene) %>%
        left_join(gpacl) %>%
        mutate(gene = factor(gene, gene_order)) %>%
        ggplot() +
        geom_segment(aes(x = gene, xend = gene, y = or-0.01, yend = or+0.01), alpha = 0.4) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
        ) +
        guides() +
        labs(x = "accessory gene", y = "log(ad/bc), if ad > bc")

    if (by_replicon) {
        return(p + facet_grid(.~replicon_type, scales = "free_x", space = "free_x"))
    } else return(p)

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
        #left_join(isolates_tax) %>%
        #mutate(contig_species = str_remove(contig_species, "E. ")) %>%
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

# 2. Elevation medicae  ----
set_name = "elev_med"
tt <- read_gpas(set_name)
n_all <- nrow(tt$gene_order) # number of all genes
n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes

# Genes that are assigned to diffeent contigs in different genomes
list_wgpa <- tt$gpacl %>%
    select(gene, genome_id, replicon_type) %>%
    distinct(gene, genome_id, .keep_all = T) %>%
    group_by(gene) %>%
    summarize(is_all_same = length(unique(replicon_type[!is.na(replicon_type)])) %in% c(0,1)) %>%
    filter(!is_all_same)

p1 <- plot_heatmap(tt$gpa, tt$gpatl, tt$gpacl, list_wgpa) + ggtitle(paste0(n_accessory, " accessory genes"))
tidied_fisher <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/tidied_fisher.csv"))
p2 <- plot_fisher(tidied_fisher, tt$gpa, tt$gpatl, tt$gpacl, list_wgpa) + ggtitle("Odds ratio")
p <- plot_grid(p1, p2, nrow = 2, axis = "lrt", align = "vh")
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-01-gpa_heatmap.png"), p, width = 10, height = 6)

p <- plot_gfs(tt$gpa)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-02-gene_frequency_spectrum.png"), p, width = 5, height = 4)
p <- plot_singleton(tt$gpa, tt$gpatl)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-03-singletons.png"), p, width = 10, height = 8)
p <- plot_bcsm(tt$sml)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-04-gpa_sm.png"), p, width = 10, height = 8)

# tt$gpacl %>%
#     select(gene, genome_id, replicon_type) %>%
#     distinct(gene, genome_id, .keep_all = T) %>%
#     pivot_wider(names_from = genome_id, values_from = replicon_type, values_fill = NA) %>%
#     filter(gene == "xerD_2~~~xerC_1") %>%
#     view



# 3. Urbanization meliloti  ----
set_name = "urbn_mel"
tt <- read_gpas(set_name)
n_all <- nrow(tt$gene_order) # number of all genes
n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes


p1 <- plot_heatmap(tt$gpa, tt$gpatl, tt$gpacl) + ggtitle(paste0(n_accessory, " accessory genes"))
tidied_fisher <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/tidied_fisher.csv"))
p2 <- plot_fisher(tidied_fisher, tt$gpa, tt$gpatl, tt$gpacl) + ggtitle("Odds ratio")
p <- plot_grid(p1, p2, nrow = 2, axis = "lrt", align = "vh")
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-01-gpa_heatmap.png"), p, width = 10, height = 6)


p <- plot_gfs(tt$gpa)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-02-gene_frequency_spectrum.png"), p, width = 5, height = 4)
p <- plot_singleton(tt$gpa, tt$gpatl)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-03-singletons.png"), p, width = 10, height = 8)
p <- plot_bcsm(tt$sml)
ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-04-gpa_sm.png"), p, width = 10, height = 8)

tidied_fisher %>%
    ungroup() %>%
    arrange(desc(or)) %>%
    slice_max(or, prop = 0.01) %>%
    filter(!str_detect(gene, "group")) %>%
    write_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/top_gene_or.csv"))


























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

