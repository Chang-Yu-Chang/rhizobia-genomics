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
make_gene_fst <- function (gene_lengths, per_gene_fst, gpacl) {
    # Use the first genome
    gene_replicon <- filter(gpacl, genome_id == gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))

    gene_lengths %>%
        filter(gene %in% per_gene_fst$gene) %>%
        left_join(per_gene_fst) %>% # those genes without Fst do not have SNPs
        left_join(gene_replicon) %>% # get the replicon where the gene is
        group_by(replicon_type) %>%
        mutate(loc_start = cumsum(sequence_length)) %>%
        mutate(loc_start = lag(loc_start)) %>%
        replace_na(list(loc_start = 0, Gst_est = 0)) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce")))

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
plot_gene_fst <- function (gene_fst, labelled_genes = "fix|nod|nif") {
    #xx_snps <- distinct(ungroup(xx_gn), n_snps)
    gene_fst %>%
        mutate(labelled_col = ifelse(str_detect(gene, labelled_genes), T, F)) %>%
        mutate(loc_start = loc_start / 10^6) %>%
        ggplot() +
        geom_segment(aes(x = loc_start, xend = loc_start, y = fst-0.01, yend = fst+0.01, color = labelled_col, alpha = labelled_col)) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3)) +
        scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black")) +
        facet_grid(.~replicon_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme() +
        guides() +
        labs(x = "core genome (Mbp)", y = "fst") +
        ggtitle(label = "", subtitle = labelled_genes)

}
plot_snps_fst <- function (snp_fst) {
    snp_fst %>%
        mutate(loc_start = loc_start / 10^6) %>%
        ggplot() +
        geom_segment(aes(x = loc_start, xend = loc_start, y = fst-0.01, yend = fst+0.01), alpha = 0.4) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        facet_grid(.~replicon_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(
            panel.grid.major =
        ) +
        guides() +
        labs(x = "core genome (Mbp)", y = "fst")
}
plot_box <- function (gene_fst, labelled_genes = "fix|nod|nif") {

}
slice_top_genes <- function (gene_fst, prop = 0.01) {
    # Genes with top ORs
    gene_fst %>%
        ungroup() %>%
        rename_with(.col = starts_with("fst."), function (x) str_remove(x, "fst.")) %>%
        arrange(desc(fst)) %>%
        slice_max(fst, prop = 0.01) %>%
        # remove unknown genes
        #filter(!str_detect(gene, "group")) %>%
        arrange(replicon_type)
}

#
set_name = "elev_med"
#set_name = "urbn_mel"
tt <- read_gpas(set_name)
ff <- read_fsts(set_name)

gene_fst <- make_gene_fst(ff$gene_lengths, ff$per_gene_fst, tt$gpacl)
snp_fst <- make_snp_fst(gene_fst, ff$per_locus_fst)

nrow(gene_fst) # number of single copy core genes
nrow(snp_fst) # number of snps

p1 <- plot_gene_fst(gene_fst, "fix|nod|nif") + ggtitle(paste0(set_name, ": ", nrow(gene_fst), " single copy core genes"))
p2 <- plot_snps_fst(snp_fst) + ggtitle(paste0(set_name, ": ", nrow(snp_fst), " SNPs"))
p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-01-fst.png"), p, width = 10, height = 8)


# Site frequence spectrum
# snp_fst %>%
#     ggplot() +
#     geom_histogram(aes(h1)) +
#     theme_bw() +
#     theme() +
#     guides() +
#     labs()


if (F) {
    labelled_genes = "fix|nod|nif|noe|fdx"

    gene_fst %>%
        filter(str_detect(gene, labelled_genes))

    tt$gpa %>%
        filter(str_detect(gene, labelled_genes)) %>%
        view
    # acce_fst %>%
    #     filter(str_detect(gene, labelled_genes))

    gene_fst %>%
        mutate(targeted = str_detect(gene, labelled_genes)) %>%
        #filter(replicon_type == "p")
        ggplot() +
        geom_boxplot(aes(x = targeted, y = fst)) +
        geom_jitter(aes(x = targeted, y = fst), shape = 21, width = 0.1) +
        facet_grid(.~replicon_type) +
        theme_bw() +
        theme() +
        guides() +
        labs()

}
