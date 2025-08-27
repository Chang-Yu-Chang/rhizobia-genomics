#' random forest importance score

library(tidyverse)
library(cowplot)
library(vegan)
library(randomForest)
source(here::here("metadata.R"))


isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>%
    mutate(population = ifelse(population == "VA", "Virginia", "Pennsylvania"))
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    mutate(population = ifelse(population == "VA", "Virginia", "Pennsylvania"))
tt <- read_gpas()

tb <- tt$gd %>%
    filter(str_detect(gene, "dna|grp|gro|rpo|clp|rec|uvr")) %>%
    mutate(ge = str_sub(gene, 1, 5) %>% str_remove("_\\d$|_$")) %>%
    mutate(g = str_sub(ge, 1, 3)) %>%
    select(g, ge, gene, genome_id) %>%
    left_join(select(iso, genome_id, contig_species))

tbm <- tb %>%
    filter(contig_species %in% c("S. meliloti", "S. medicae")) %>%
    group_by(contig_species, genome_id, ge) %>%
    count() %>%
    pivot_wider(names_from = ge, values_from = n, values_fill = 0)
dm <- vegdist(tbm[,-c(1,2)], method = "bray")

# Panel A. Which genes differ ---
rf <- randomForest(
    x = tbm[,-c(1,2)],
    y = as.factor(tbm$contig_species),
    importance = TRUE,
    ntree = 500
)

p1 <- importance(rf) %>%
    as_tibble() %>%
    mutate(ge = rownames(importance(rf))) %>%
    mutate(percentage = MeanDecreaseGini / sum(MeanDecreaseGini)* 100) %>%
    arrange(MeanDecreaseGini) %>%
    mutate(ge = factor(ge, ge)) %>%
    ggplot() +
    geom_col(aes(x = ge, y = percentage)) +
    coord_flip(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs(x = "", y = "Scaled relative importance (%)")


# Panel B. enriched genes in what direction??
set.seed(123)
ord <- capscale(dm ~ 1)   # unconstrained ordination (PCoA)
genes_focus <- c("recA", "uvrD", "dnaE2")
ef <- envfit(ord, tbm[, genes_focus], permutations = 999)

# site scores
scores_sites <- scores(ord)$sites %>%
    as_tibble() %>%
    mutate(genome_id = tbm$genome_id) %>%
    left_join(select(tbm, genome_id, contig_species))

# gene loadings
ef_df <- as.data.frame(ef$vectors$arrows * sqrt(ef$vectors$r))
ef_df$gene <- rownames(ef_df)

p2 <- scores_sites %>%
    ggplot() +
    geom_point(aes(x = MDS1, y = MDS2, color = contig_species), size = 3) +
    stat_ellipse(aes(x = MDS1, y = MDS2, color = contig_species), level = 0.95) +
    geom_segment(data = ef_df, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), arrow = arrow(length = unit(0.2,"cm")), color = "black") +
    geom_text(data = ef_df, aes(x = MDS1, y = MDS2, label = gene), size = 4, vjust = -0.5) +
    scale_color_manual(values = species_colors) +
    coord_fixed(clip = "off") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.line = element_line(color = 1),
        legend.title = element_blank(),
        legend.position = "top"
    )

# ----
p <- plot_grid(p1, plot_grid(p2, NULL, rel_heights = c(1, .3), ncol = 1), labels = c("A", "B"), scale = .95, rel_widths = c(1, 1.5)) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/FigS6.png"), p, width = 8, height = 6)
