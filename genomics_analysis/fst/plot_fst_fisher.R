#' This script plot the Fst per gene

renv::load()
library(tidyverse)
library(cowplot)
library(ggsci)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
set1_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/set1_fst.csv"))
set1_pop_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/set1_pop_fst.csv"))
set1_fst_sccg <- set1_fst %>%
    filter(metric == "Gst_est", singlecopy) %>%
    select(-metric, -singlecopy) %>%
    pivot_longer(-gene, names_to = "gradient") %>%
    mutate(gradient = ifelse(gradient == "elev", "elevation", "urbanization"))

# 1. Plot the Fst within gradients ----
n_sccg <- distinct(set1_fst, gene, singlecopy) %>% pull(singlecopy) %>% sum

p <- set1_fst_sccg %>%
    ggplot() +
    geom_histogram(aes(x = value, fill = gradient), alpha = .5, position = "identity", binwidth = 0.01) +
    scale_fill_aaas() +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(.8,.8),
        legend.background = element_rect(color = "grey10", fill = "white")
    ) +
    guides() +
    labs(x = "Nei's G_st", title = paste0(n_sccg, " single-copy core genes"))
ggsave(paste0(folder_data, "genomics_analysis/fst/01-fst_sccg.png"), p, width = 4, height = 4)


# 2. Plot the Fst between all populations across the two gradients ----
n_sccg <- distinct(set1_fst, gene, singlecopy) %>% pull(singlecopy) %>% sum
p <- set1_pop_fst %>%
    filter(singlecopy) %>%
    select(-singlecopy) %>%
    rename(gradient1 = pop1, gradient2 = pop2) %>%
    #pivot_longer(cols = -c(gene, Gst_est), names_to = "gradient_group", values_to = "pop_group") %>%
    ggplot() +
    geom_histogram(aes(x = Gst_est), alpha = .5, position = "identity", binwidth = 0.01) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    scale_fill_aaas() +
    coord_cartesian(clip = "off") +
    facet_grid(gradient1 ~ gradient2) +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(.8,.8),
        legend.background = element_rect(color = "grey10", fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
    ) +
    guides() +
    labs(x = "Nei's G_st", title = paste0(n_sccg, " single-copy core genes"))
ggsave(paste0(folder_data, "genomics_analysis/fst/02-fst_pop_sccg.png"), p, width = 8, height = 8)


# Top genes ----
set1_fst %>%
    filter(metric == "Gst_est", singlecopy) %>%
    select(-metric, -singlecopy) %>%
    pivot_longer(-gene, names_to = "gradient") %>%
    mutate(gradient = ifelse(gradient == "elev", "elevation", "urbanization")) %>%
    group_by(gradient) %>%
    arrange(desc(value)) %>%
    slice(1:ceiling(n()/100))

# 3. GCV fisher exact historgram ----
gpa_fis <- read_csv(paste0(folder_data, "genomics_analysis/fst/gpa_fis.csv"))

p <- gpa_fis %>%
    mutate(gene = factor(gene, gene_order$gene)) %>%
    #filter(p.value < 0.05) %>%
    #filter(!gene %in% set1_fst$gene) %>%
    ggplot() +
    geom_histogram(aes(x = estimate, fill = gradient), alpha = .5, position = "identity", binwidth = 0.5) +
    scale_fill_aaas() +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(.8,.8),
        legend.background = element_rect(color = "grey10", fill = "white")
    ) +
    guides() +
    labs(x = "Fisher's exact", title = paste0(nrow(gene_order), " genes in pangenome"))

ggsave(paste0(folder_data, "genomics_analysis/fst/03-fis_gcv.png"), p, width = 4, height = 4)

# 4. Fst manhattan plot ----
gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv")) %>%
    mutate(gene_id = factor(1:n()))
gpa_fis <- read_csv(paste0(folder_data, "genomics_analysis/fst/gpa_fis.csv"))
grt <- gpacl %>%
    drop_na(replicon_type) %>%
    distinct(gene, replicon_type) %>%
    distinct(gene, .keep_all = T)

tb1 <- set1_fst_sccg %>%
    left_join(grt) %>%
    left_join(gene_order) %>%
    filter(gradient == "elevation") %>%
    mutate(gene = factor(gene, gene_order$gene))
tb2 <- gpa_fis %>%
    select(gene, gradient, estimate, p.value) %>%
    left_join(grt) %>%
    left_join(gene_order) %>%
    filter(gradient == "elevation") %>%
    mutate(gene = factor(gene, gene_order$gene))

# gpa
background_df <- tibble(site_group = c("high elevation", "low elevation", "suburban", "urban"))

p1 <- gpacl %>%
    #filter(genome_id != "g42") %>%
    #drop_na(replicon_type) %>%
    left_join(isolates) %>%
    filter(population == "VA") %>%
    mutate(gene = factor(gene, gene_order$gene)) %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    left_join(gene_order) %>%
    ggplot() +
    geom_tile(aes(x = gene_id, y = genome_id, fill = site_group)) +
    scale_fill_manual(values = alpha(site_group_colors[1:2], 0.8))+
    #scale_x_discrete(expand = c(0,0), drop = F) +
    #scale_x_continuous(limits = c(1, nrow(gene_order)), breaks = c(1, seq(1000, 26000, 1000)), expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0), position = "right") +
    facet_grid(site_group ~ replicon_type, scales = "free", space = "free", switch = "y") +
    #facet_grid(site_group ~., scales = "free", space = "free", switch = "y" ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        #axis.title.x = element_text(size = 15),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.spacing.y = unit(0, "mm"),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0, size = 10),
        strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0),
        strip.text.y.left = element_text(angle = 75, hjust = .8, vjust = .5, size = 10),
        strip.clip = "off",
        plot.background = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "gene", y = "genome")

## Fst

p2 <- tb1 %>%
    mutate(istop1 = value > quantile(value, 0.99)) %>%
    #drop_na(replicon_type) %>%
    ggplot() +
    geom_point(aes(x = gene_id, y = value, color = istop1), size = .5) +
    #scale_x_continuous(limits = c(1, nrow(gene_order)), breaks = c(1, seq(1000, 26000, 1000)), expand = c(0,0)) +
    facet_grid(~replicon_type, scales = "free", space = "free", switch = "y") +
    scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "red")) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0, size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    guides() +
    labs(x = "gene", y = "Nei's G_st", title = "777 single-copy core genes")

## Fisher manhattan
p3 <- tb2 %>%
    mutate(signif = p.value < 0.05) %>%
    #drop_na(replicon_type) %>%
    ggplot() +
    geom_point(aes(x = gene_id, y = estimate, color = signif), size = .5) +
    scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "red")) +
    facet_grid(~replicon_type, scales = "free", space = "free", switch = "y") +
    #scale_x_continuous(breaks = c(1, seq(1000, 26000, 1000)), expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0, size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    guides() +
    labs(x = "gene", y = "Fisher odds ratio", title = "Gene presence/absence")

## -log P
p4 <- tb2 %>%
    mutate(signif = p.value < 0.05) %>%
    #drop_na(replicon_type) %>%
    ggplot() +
    geom_point(aes(x = gene_id, y = -log(p.value), color = signif), size = .5) +
    geom_hline(yintercept = -log(0.05), linetype = 2, color = "maroon") +
    scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "red")) +
    facet_grid(~replicon_type, scales = "free", space = "free", switch = "y") +
    #scale_x_continuous(breaks = c(1, seq(1000, 26000, 1000)), expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0, size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    guides() +
    labs(x = "gene", y = "-log(P)", title = "Gene presence/absence")

p <- plot_grid(p1,p2,p3,p4, ncol = 1, axis = "lrt", align = "vh", labels = LETTERS[1:4]) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(paste0(folder_data, "genomics_analysis/fst/04-elevation.png"), p, width = 10, height = 10)

if (F) {




isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv"))
gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpacl.csv")) %>%
    separate(contig_id, into = c("genome_id", "temp1", "temp2"), remove = F) %>%
    select(-temp1, -temp2)
background_df <- tibble(site_group = c("high elevation", "low elevation", "suburban", "urban"))

p_gpam <- gpacl %>%
    filter(genome_id != "g42") %>%
    drop_na(replicon_type) %>%
    left_join(isolates) %>%
    mutate(gene = factor(gene, gene_order$gene)) %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = site_group), alpha = 0.2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
    scale_fill_manual(values = site_group_colors) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0), position = "right") +
    facet_grid(site_group ~ replicon_type, scales = "free", space = "free", switch = "y" ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        strip.background = element_blank(),
        panel.spacing.y = unit(c(0,3,0), "mm"),
        strip.text.x = element_text(angle = 15, hjust = 0, vjust = 0, size = 15),
        #strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0),
        strip.text.y.left = element_text(angle = 75, hjust = .8, vjust = .5, size = 15),
        strip.clip = "off",
        plot.background = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome")

nrow(gene_order) # 26504 genes in the pangenome
}
