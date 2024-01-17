#' This script joins the phenotypic distance

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    library(ggsci)
    source(here::here("analysis/00-metadata.R"))
})


dists_genetic <- read_csv(paste0(folder_data, 'temp/19-dists_genetic.csv'), show_col_types = F)
dists_trait <- read_csv(paste0(folder_data, "temp/29-dists_trait.csv"), show_col_types = F)

# 0. ----
dists <- dists_genetic %>%
    left_join(dists_trait) %>%
    select(-row_name, -col_name) %>%
    rename_with(~str_replace(.x, "distance", "d"), starts_with("distance"))


write_csv(dists, paste0(folder_data, "temp/31-dists.csv"))

# 1.  one for example ----
p <- dists %>%
    ggplot() +
    geom_point(aes(x = d_kmer, y = d_growth), shape = 21) +
    geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/31-01-kmer_vs_growth.png"), p, width = 4, height = 4)

# 2. boxplot ----
dists_long <- dists %>%
    select(-row_name, -col_name) %>%
    pivot_longer(cols = c(-genome_id1, -genome_id2), names_to = "d_type", names_prefix = "d_")

p <- dists_long %>%
    drop_na(value) %>%
    mutate(d_type = factor(d_type, c("kmer", "jaccard", "fluidity", "growth", "symbiosis"))) %>%
    ggplot() +
    geom_boxplot(aes(x = d_type, y = value)) +
    geom_jitter(aes(x = d_type, y = value), width = 0.3, shape = 21) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/31-02-boxplot.png"), p, width = 6, height = 4)


# 3. composite distance comparison ----
plot_dots <- function (tb, var1, var2) {
    tb %>%
        ggplot() +
        geom_point(aes(x = {{var1}}, y = {{var2}}), shape = 21) +
        geom_abline(intercept = 0, slope = 1, color = "maroon", linetype = 2) +
        scale_x_continuous(limits = c(0,1)) +
        scale_y_continuous(limits = c(0,1)) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        guides() +
        labs()
}

p1 <- plot_dots(dists, d_kmer, d_growth)
p2 <- plot_dots(dists, d_jaccard, d_growth)
p3 <- plot_dots(dists, d_fluidity, d_growth)
p4 <- plot_dots(dists, d_kmer, d_symbiosis)
p5 <- plot_dots(dists, d_jaccard, d_symbiosis)
p6 <- plot_dots(dists, d_fluidity, d_symbiosis)


p <- plot_grid(p1, p2, p3, p4, p5, p6, byrow = T)
ggsave(paste0(folder_data, "temp/31-03-pairwise_distance.png"), p, width = 10, height = 6)

# Stat
cor.test(dists$d_kmer, dists$d_growth)
cor.test(dists$d_jaccard, dists$d_growth)
cor.test(dists$d_fluidity, dists$d_growth)
cor.test(dists$d_kmer, dists$d_symbiosis)
cor.test(dists$d_jaccard, dists$d_symbiosis)
cor.test(dists$d_fluidity, dists$d_symbiosis)
# Basically all significant










