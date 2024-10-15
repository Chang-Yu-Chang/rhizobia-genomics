#' This script plot the Fst per gene and per snp

library(tidyverse)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
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
read_fsts <- function (set_name) {
    gene_wide_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_wide_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))
    return(list(gene_wide_fst = gene_wide_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths))
}
make_gene_fst <- function (gene_lengths, gene_wide_fst, gpacl) {
    # Use the first genome
    gene_replicon <- filter(gpacl, genome_id == gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))
    gene_lengths %>%
        filter(gene %in% gene_wide_fst$gene) %>%
        left_join(gene_wide_fst) %>% # those genes without Fst do not have SNPs
        left_join(gene_replicon) %>% # get the replicon where the gene is
        group_by(replicon_type) %>%
        mutate(loc_start = cumsum(sequence_length)) %>%
        mutate(loc_start = lag(loc_start)) %>%
        replace_na(list(loc_start = 0, Gst_est = 0)) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "psymA like", "psymB like", "others")))

}
make_snp_fst <- function (gene_fst, per_locus_fst) {
    gene_fst %>%
        select(gene, replicon_type, loc_start) %>%
        group_by(replicon_type) %>%
        right_join(per_locus_fst) %>%
        ungroup() %>%
        arrange(replicon_type, loc_start) %>%
        mutate(locus_id = 1:n()) %>%
        mutate(loc_start = loc_start + location)
}
plot_gene_fst <- function (gene_fst) {
    gene_fst %>%
        mutate(loc_start = loc_start / 10^6) %>%
        ggplot() +
        geom_segment(aes(x = loc_start, xend = loc_start, y = fst-0.01, yend = fst+0.01), alpha = 0.4) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        scale_y_continuous(limits = c(-0.01, 1.01)) +
        facet_grid(.~replicon_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme() +
        guides() +
        labs(x = "core genome (Mbp)", y = "fst")

}
plot_snps_fst <- function (snp_fst) {
    snp_fst %>%
        mutate(loc_start = loc_start / 10^6) %>%
        ggplot() +
        geom_segment(aes(x = loc_start, xend = loc_start, y = fst-0.01, yend = fst+0.01), alpha = 0.4) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        scale_y_continuous(limits = c(-0.01, 1.01)) +
        facet_grid(.~replicon_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(
            panel.grid.major =
        ) +
        guides() +
        labs(x = "core genome (Mbp)", y = "fst")
}
slice_top_genes <- function (gene_fst, prop = 0.01) {
    # Genes with top ORs
    gene_fst %>%
        ungroup() %>%
        arrange(desc(fst)) %>%
        slice_max(fst, prop = 0.01) %>%
        # remove unknown genes
        #filter(!str_detect(gene, "group")) %>%
        arrange(replicon_type)
}

# Elevation medicae  ----
set_name = "elev_med"
tt <- read_gpas(set_name)
ff <- read_fsts(set_name)

gene_fst <- make_gene_fst(ff$gene_lengths, ff$gene_wide_fst, tt$gpacl)
snp_fst <- make_snp_fst(gene_fst, ff$per_locus_fst)
top_gene_fst <- slice_top_genes(gene_fst)
write_csv(top_gene_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/top_gene_fst.csv"))

nrow(gene_fst) # number of single copy core genes
nrow(snp_fst) # number of snps

p1 <- plot_gene_fst(gene_fst) + ggtitle(paste0(set_name, ": ", nrow(gene_fst), " single copy core genes"))
p2 <- plot_snps_fst(snp_fst) + ggtitle(paste0(set_name, ": ", nrow(snp_fst), " SNPs"))
p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-01-fst.png"), p, width = 10, height = 8)




# Urbanization meliloti  ----
set_name = "urbn_mel"
tt <- read_gpas(set_name)
ff <- read_fsts(set_name)

gene_fst <- make_gene_fst(ff$gene_lengths, ff$gene_wide_fst, tt$gpacl)
snp_fst <- make_snp_fst(gene_fst, ff$per_locus_fst)
top_gene_fst <- slice_top_genes(gene_fst)
write_csv(top_gene_fst, paste0(folder_data, "genomics_analysis/fst/", set_name, "/top_gene_fst.csv"))

nrow(gene_fst) # number of single copy core genes
nrow(snp_fst) # number of snps

p1 <- plot_gene_fst(gene_fst) + ggtitle(paste0(set_name, ": ", nrow(gene_fst), " single copy core genes"))
p2 <- plot_snps_fst(snp_fst) + ggtitle(paste0(set_name, ": ", nrow(snp_fst), " SNPs"))
p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-01-fst.png"), p, width = 10, height = 8)


nspg <- snp_fst %>%
    group_by(gene, replicon_type) %>%
    summarize(n_snp_per_gene = n())


gene_fst %>%
    group_by(gene, replicon_type) %>%
    left_join(select(nspg, gene, n_snp_per_gene)) %>%
    #filter(n_snp_per_gene >= 5) %>%
    #filter(str_detect(gene, "nif|fix|nod")) %>%
    #filter(str_detect(gene, "hna")) %>%
    filter(!str_detect(gene, "group")) %>%
    filter(replicon_type == "psymB like") %>%
    arrange(desc(fst))


tt$gpa %>%
    filter(str_detect(gene, "nifH"))

pat <- "glxR"
tt$gpacl %>%
    filter(str_detect(gene, pat) | str_detect(annotation, pat)) %>%
    select(genome_id, replicon_type, gene, annotation) %>%
    distinct(replicon_type, gene, annotation)

gene_fst %>%
    filter(gene == "abo_4")

    # ggplot() +
    # geom_histogram(aes(x = n_snp_per_gene), binwidth = 1) +
    # theme_bw() +
    # theme() +
    # guides() +
    # labs()



