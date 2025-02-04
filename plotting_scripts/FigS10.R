#' This script plot the Fst per gene and per snp

library(tidyverse)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
read_gcv_fsts <- function (set_name) {
    per_acce_fst <- read_csv(paste0(folder_data, "genomics_analysis/gcv_fst/", set_name,"/per_acce_fst.csv"))
    per_genome_fst <- read_csv(paste0(folder_data, "genomics_analysis/gcv_fst/", set_name,"/per_genome_fst.csv"))
    return(list(per_acce_fst = per_acce_fst, per_genome_fst = per_genome_fst))
}

plot_snps_fst_hist <- function (ff) {
    ff$per_acce_fst %>%
        drop_na(Gprime_st) %>%
        #arrange(desc(Gprime_st))
        ggplot() +
        geom_histogram(aes(x = Gprime_st), binwidth = .05, color = 1, fill = "grey90") +
        scale_x_continuous(limits = c(-.4, 1.05), breaks = seq(-1, 1, .5)) +
        scale_y_continuous(expand = c(.02,0), limits = c(0, 2000)) +
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
plot_snps_fst <- function (tt, ff) {
    gene_replicon <- filter(tt$gpacl, genome_id == tt$gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others")) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce", "others")))
    ff$per_acce_fst %>%
        left_join(gene_replicon) %>%
        replace_na(list(replicon_type = "others")) %>%
        drop_na(Gprime_st) %>%
        ggplot() +
        geom_point(aes(x = gene, y = Gprime_st), shape = 21, size = .5) +
        scale_y_continuous(limits = c(-.4, 1.05), breaks = seq(-1, 1, .5)) +
        facet_grid2(~replicon_type, scales = "free_x", space = "free_x", strip = strip_vanilla(clip = "off")) +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            strip.text.x = element_text(hjust = 0, vjust = 0, angle = 45),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides() +
        labs(x = "gene cluster", y = "G'st")
}
set_name <- "elev_med"
tt1 <- read_gpas(set_name)
ff1 <- read_gcv_fsts(set_name)
p_hist1 <- plot_snps_fst_hist(ff1) + ggtitle("Elevation") + theme(axis.title.x = element_blank())
p_fst1 <- plot_snps_fst(tt1, ff1)

set_name <- "urbn_mel"
tt2 <- read_gpas(set_name)
ff2 <- read_gcv_fsts(set_name)
p_hist2 <- plot_snps_fst_hist(ff2) + ggtitle("Urbanization")
p_fst2 <- plot_snps_fst(tt2, ff2)

p <- plot_grid(
    p_hist1, p_fst1, p_hist2, p_fst2,
    ncol = 2, align = "h", axis = "tb",
    rel_widths = c(1, 3), labels = LETTERS[1:4]
)

ggsave(here::here("plots/FigS10.png"), p, width = 8, height = 6)

# T test
t.test(ff1$per_acce_fst$Gprime_st, ff2$per_acce_fst$Gprime_st) # t = 11.14, df = 8317.5, p-value < 2.2e-16
