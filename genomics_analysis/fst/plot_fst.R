#' This script plot the Fst per gene and per snp

library(tidyverse)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
read_fsts <- function (set_name) {
    per_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_gene_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))
    return(list(per_gene_fst = per_gene_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths))
}
make_gene_fst <- function (ff, tt) {
    # Use the first genome
    gene_replicon <- filter(tt$gpacl, genome_id == tt$gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))

    ff$gene_lengths %>%
        filter(gene %in% ff$per_gene_fst$gene) %>%
        left_join(ff$per_gene_fst) %>% # those genes without Fst do not have SNPs
        left_join(gene_replicon) %>% # get the replicon where the gene is
        group_by(replicon_type) %>%
        mutate(loc_start = cumsum(sequence_length)) %>%
        mutate(loc_start = lag(loc_start)) %>%
        replace_na(list(loc_start = 0, Gst_est = 0)) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        ungroup()

}

make_gene_fst_long <- function (gene_fst) {
    gene_fst %>%
        pivot_longer(cols = c(Hs, Ht, Gst_est, Gprime_st, D_het, D_mean), names_to = "Fst") %>%
        mutate(Fst = factor(Fst, c("Hs", "Ht", "Gst_est", "Gprime_st", "D_het", "D_mean")))
}
make_snp_fst <- function (gene_fst, ff) {
    gene_fst %>%
        select(gene, replicon_type, loc_start) %>%
        group_by(replicon_type) %>%
        right_join(ff$per_locus_fst) %>%
        ungroup() %>%
        arrange(replicon_type, loc_start) %>%
        mutate(locus_id = 1:n()) %>%
        mutate(loc_start = loc_start + location)
}
make_snp_fst_long <- function (snp_fst) {
    snp_fst %>%
        pivot_longer(cols = c(Hs, Ht, Gst, Gprime_st, D), names_to = "Fst") %>%
        mutate(Fst = factor(Fst, c("Hs", "Ht", "Gst", "Gprime_st", "D")))
}
plot_gene_fst <- function (gene_fst_long) {
    gene_fst_long %>%
        mutate(loc_start = loc_start / 10^6) %>%
        drop_na(value) %>%
        filter(replicon_type != "pAcce") %>%
        ggplot() +
        geom_point(aes(x = loc_start, y = value), size = .5, shape = 21) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        facet_grid(Fst~replicon_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme() +
        guides() +
        labs(x = "core genome (Mbp)") +
        ggtitle(label = "")

}
plot_snps_fst <- function (snp_fst_long) {
    snp_fst_long %>%
        mutate(loc_start = loc_start / 10^6) %>%
        ggplot() +
        geom_point(aes(x = loc_start, y = value), shape = 21, size = .5) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        facet_grid(Fst~replicon_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme() +
        guides() +
        labs(x = "core genome (Mbp)", y = "fst")
}
plot_box <- function (gene_fst, labelled_genes = "fix|nod|nif") {

}
plot_fst_wrapper <- function (set_name) {
    #set_name = "elev_med"
    #set_name = "urbn_mel"
    tt <- read_gpas(set_name)
    ff <- read_fsts(set_name)

    gene_fst <- make_gene_fst(ff, tt)
    gene_fst_long <- make_gene_fst_long(gene_fst)
    snp_fst <- make_snp_fst(gene_fst, ff)
    snp_fst_long <- make_snp_fst_long(snp_fst)
    nrow(gene_fst) # number of single copy core genes
    nrow(snp_fst) # number of snps in these genes

    p <- plot_gene_fst(gene_fst_long)
    ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-01-gene_fst.png"), p, width = 10, height = 8)

    gene_fst_long %>%
        group_by(replicon_type, Fst) %>%
        arrange(desc(value)) %>%
        slice_max(value, prop = .01) %>%
        filter(replicon_type == "pSymA")

    p <- plot_snps_fst(snp_fst_long)
    ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-02-snp_fst.png"), p, width = 10, height = 8)

    snp_fst_long %>%
        group_by(replicon_type, Fst) %>%
        arrange(desc(value)) %>%
        slice_max(value, prop = .01) %>%
        filter(!Fst %in% c("Hs", "Ht"))
        #filter(replicon_type == "pSymA") %>%



    #p1 <- plot_gene_fst(gene_fst, "fix|nod|nif") + ggtitle(paste0(set_name, ": ", nrow(gene_fst), " single copy core genes"))
    #p2 <- plot_snps_fst(snp_fst) + ggtitle(paste0(set_name, ": ", nrow(snp_fst), " SNPs"))
    #p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
    #ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-01-fst.png"), p, width = 10, height = 8)


}
plot_fst_wrapper("elev_med")
plot_fst_wrapper("urbn_mel")
