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
plot_snps_fst <- function (snp_fst_long) {
    snp_fst_long %>%
        mutate(loc_start = loc_start / 10^6) %>%
        filter(replicon_type != "pAcce") %>%
        filter(Fst == "Gprime_st") %>%
        ggplot() +
        geom_point(aes(x = loc_start, y = value), shape = 21, size = .5) +
        scale_x_continuous(breaks = seq(0, 4, 0.5)) +
        scale_y_continuous(limits = c(-.4, 1), breaks = seq(-1, 1, .5)) +
        facet_grid2(~replicon_type, scales = "free_x", space = "free_x") +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            panel.grid.major = element_line(color = "grey90")
        ) +
        guides() +
        labs(x = "core genome (Mbp)", y = "G'st")
}
plot_snps_fst_hist <- function (snp_fst_long) {
    # Genome
    snp_fst_long %>%
        filter(Fst == "Gprime_st") %>%
        drop_na(value) %>%
        ggplot() +
        geom_histogram(aes(x = value), binwidth = .05, color = 1, fill = "grey90") +
        scale_x_continuous(limits = c(-.4, 1), breaks = seq(-1, 1, .5)) +
        scale_y_continuous(expand = c(.02,0), limits = c(0, 6000)) +
        theme_classic() +
        coord_flip(clip = "off") +
        theme(
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_line(color = "grey90")
        ) +
        guides() +
        labs(x = "G'st")
}


plot_fst <- function (set_name) {
    tt <- read_gpas(set_name)
    ff <- read_fsts(set_name)
    gene_fst <- make_gene_fst(ff, tt)
    snp_fst <- make_snp_fst(gene_fst, ff)
    snp_fst_long <- make_snp_fst_long(snp_fst)
    nrow(gene_fst) # number of single copy core genes with SNPs
    nrow(snp_fst) # number of snps in these genes

    gra <- ifelse(set_name == "elev_med", "Elevation", "Urbanization")

    p_hist <- plot_snps_fst_hist(snp_fst_long) + ggtitle(gra)
    p_fst <- plot_snps_fst(snp_fst_long)
    return(list(p_hist, p_fst))
}
p1 <- plot_fst("elev_med")
p2 <- plot_fst("urbn_mel")

p <- plot_grid(
    p1[[1]], p1[[2]], p2[[1]], p2[[2]],
    ncol = 2, align = "h", axis = "tb",
    rel_widths = c(1, 3), labels = LETTERS[1:4]
)

ggsave(here::here("plots/FigS8.png"), p, width = 8, height = 6)

# T test
get_snp_get <- function (set_name) {
    tt <- read_gpas(set_name)
    ff <- read_fsts(set_name)
    gene_fst <- make_gene_fst(ff, tt)
    snp_fst <- make_snp_fst(gene_fst, ff)
    return(snp_fst)
}

snp_fst1 <- get_snp_get("elev_med")
snp_fst2 <- get_snp_get("urbn_mel")

t.test(snp_fst1$Gprime_st, snp_fst2$Gprime_st)

# Numbers
nrow(snp_fst1) # 26042 snps
length(unique(snp_fst1$gene)) # 3169 single-copy core genes with snps
tt <- read_gpas("elev_med")
nrow(tt$list_sccg) # 4702

nrow(snp_fst2) # 41060 snps
length(unique(snp_fst2$gene)) # 3609 single-copy core genes with snps
tt <- read_gpas("urbn_mel")
nrow(tt$list_sccg) # 4417

