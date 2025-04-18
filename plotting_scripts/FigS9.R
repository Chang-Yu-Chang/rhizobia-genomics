#' GO terms for SNPs

library(tidyverse)
library(cowplot)
source(here::here("metadata.R"))

read_goenrich <- function (set_name) {
    goenrich_bygene <- read_csv(paste0(folder_data, "genomics_analysis/go/", set_name, "/goenrich_bygene.csv")) %>%
        mutate(goidterm = paste0(go_id, " ", term))
    goenrich_bysnp <- read_csv(paste0(folder_data, "genomics_analysis/go/", set_name, "/goenrich_bysnp.csv")) %>%
        mutate(goidterm = paste0(go_id, " ", term))
    top_genes_bygene <- read_csv(paste0(folder_data, "genomics_analysis/go/", set_name, "/top_genes_bygene.csv"))
    top_genes_bysnp <- read_csv(paste0(folder_data, "genomics_analysis/go/", set_name, "/top_genes_bysnp.csv"))

    return(list(goenrich_bygene = goenrich_bygene, goenrich_bysnp = goenrich_bysnp, top_genes_bygene = top_genes_bygene, top_genes_bysnp = top_genes_bysnp))
}
plot_go_num <- function (goenrich) {
    goenrich %>%
        mutate(goidterm = factor(goidterm, rev(goenrich$goidterm))) %>%
        ggplot() +
        geom_col(aes(x = goidterm, y = annotated, fill = category)) +
        coord_flip() +
        #scale_y_continuous(limits = c(0, 75)) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.y = element_blank()
        ) +
        guides() +
        labs(x = "", title = set_name)
}
plot_go_p <- function (goenrich) {
    goenrich %>%
        mutate(goidterm = factor(goidterm, rev(goenrich$goidterm))) %>%
        ggplot() +
        geom_col(aes(x = goidterm, y = -log(as.numeric(classic_fisher), 10), fill = category)) +
        geom_hline(yintercept = -log(0.05, 10)) +
        coord_flip() +
        scale_x_discrete(position = "top") +
        theme_bw() +
        theme(
            axis.text.y = element_text(hjust = 0),
            legend.position = "top"
        ) +
        guides() +
        labs(x = "", y = "-log(p)")
}

set_name = "elev_med"
ee1 <- read_goenrich(set_name)
p1a <- plot_go_num(ee1$goenrich_bygene) + ggtitle("Elevation")
p1b <- plot_go_p(ee1$goenrich_bygene)

set_name = "urbn_mel"
ee2 <- read_goenrich(set_name)
p2a <- plot_go_num(ee2$goenrich_bygene) + ggtitle("Urbanization")
p2b <- plot_go_p(ee2$goenrich_bygene) + guides(fill = "none")

p <- plot_grid(
    p1a, p1b, p2a, p2b,
    nrow = 2, align = "h", axis = "t", rel_widths = c(.3, .8),
    labels = c("A", "", "B", "")
)

ggsave(here::here("plots/FigS9.png"), p, width = 8, height = 10)

# Top genes
ee1$top_genes_bygene %>% arrange(desc(Gprime_st)) # yhdY
nrow(ee1$top_genes_bygene) # 1872 core genes with GO terms
sum(ee1$top_genes_bygene$tops) # 103 top genes

ee2$top_genes_bygene %>% arrange(desc(Gprime_st)) # paoC and nanR
nrow(ee2$top_genes_bygene) # 2200 core genes with GO terms
sum(ee2$top_genes_bygene$tops) # 117 top genes
