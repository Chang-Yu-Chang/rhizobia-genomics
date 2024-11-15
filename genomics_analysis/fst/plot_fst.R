#' This script plot the Fst per gene and per snp

library(tidyverse)
library(cowplot)
library(ggh4x)
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
        filter(replicon_type != "pAcce") %>%
        ggplot() +
        geom_point(aes(x = loc_start, y = value), shape = 21, size = .5) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        facet_grid(Fst~replicon_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme() +
        guides() +
        labs(x = "core genome (Mbp)", y = "fst")
}
plot_snps_fst_hist <- function (snp_fst_long) {
    # Genome
    p1 <- snp_fst_long %>%
        filter(Fst == "Gprime_st") %>%
        drop_na(value) %>%
        #pull(value) %>% range()
        ggplot() +
        geom_histogram(aes(x = value), binwidth = .1, color = 1, fill = "grey90") +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
        scale_y_continuous(expand = c(.02,0), limits = c(0, 6000)) +
        #facet_grid2(~replicon_type, axes = "all", remove_labels = T) +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        guides() +
        labs(x = "G'st")

    # Replicon
    p2 <- snp_fst_long %>%
        filter(replicon_type != "pAcce") %>%
        filter(Fst == "Gprime_st") %>%
        ggplot() +
        geom_histogram(aes(x = value), binwidth = .1, color = 1, fill = "grey90") +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
        scale_y_continuous(expand = c(.02,0), limits = c(0, 6000)) +
        facet_grid2(~replicon_type, axes = "all", remove_labels = T) +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        guides() +
        labs(x = "G'st")

    p <- plot_grid(p1, p2, rel_widths = c(1, 3), align = "h", axis = "tb")
    return(p)
}
plot_gene_fst_top <- function (gene_fst_long) {
    gene_fst_line <- gene_fst_long %>%
        filter(replicon_type != "pAcce") %>%
        group_by(replicon_type, Fst) %>%
        arrange(desc(value)) %>%
        slice_max(value, prop = .01) %>%
        filter(Fst %in% "Gprime_st") %>%
        summarize(min_value = min(value))

    gene_fst_long %>%
        mutate(loc_start = loc_start / 10^6) %>%
        filter(replicon_type != "pAcce") %>%
        filter(Fst == "Gprime_st") %>%
        ggplot() +
        geom_point(aes(x = loc_start, y = value), shape = 21, size = .5) +
        geom_hline(data = gene_fst_line, aes(yintercept = min_value), color = "red", linetype = 2) +
        scale_x_continuous(breaks = seq(0, 3, .5)) +
        scale_y_continuous(expand = c(.02,0), breaks = seq(0, 1, .5)) +
        facet_grid2(~replicon_type, axes = "all", scales = "free_x", space = "free_x") +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        guides() +
        labs(x = "core genome (Mbp)", y = "G'st")
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

    p <- plot_snps_fst(snp_fst_long)
    ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-02-snp_fst.png"), p, width = 10, height = 8)

    p <- plot_snps_fst_hist(snp_fst_long)
    ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-03-gene_fst_hist.png"), p, width = 10, height = 3)

    # Genes with top Fst
    p <- plot_gene_fst_top(gene_fst_long)
    ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-04-gene_fst_gprime.png"), p, width = 8, height = 3)

}
plot_fst_wrapper("elev_med")
plot_fst_wrapper("urbn_mel")

#
# gene_fst_long %>%
#     filter(replicon_type != "pAcce") %>%
#     mutate(loc_start = loc_start / 10^6) %>%
#     group_by(replicon_type, Fst) %>%
#     arrange(desc(value)) %>%
#     slice_max(value, prop = .01) %>%
#     filter(Fst %in% "Gprime_st") %>%
#     view
# #filter(!str_detect(gene, "group"))
#


# gene_fst_long %>%
#     group_by(replicon_type, Fst) %>%
#     arrange(desc(value)) %>%
#     slice_max(value, prop = .01) %>%
#     filter(replicon_type == "pSymA") %>%
#     view

