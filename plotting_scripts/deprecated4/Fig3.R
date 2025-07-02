#' Isolation by distance

library(tidyverse)
library(cowplot)
deprlibrary(ggh4x)
library(tidytree)
library(ggtree)
library(vegan) # for mantel test
library(lme4) # for lmer
library(car) # for anova
library(broom.mixed)
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- select(isolates, -exp_id, -genome_name)
sites_dist <- read_csv(paste0(folder_data, "phenotypes/sites/sites_dist.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
tb2 <- read_csv(paste0(folder_data, "phylogenomics_analysis/tree_distance/tb2.csv"))


# Panel A. core gene ----
nodes_to_scale <- c(38, 40, 1, 2, 41, 42, 54)
tr <- tbtr$tr[[1]]
tr <- root(tr, outgroup = "g2")
edges_to_scale <- which(tr$edge[,2] %in% nodes_to_scale)
tr$edge.length[edges_to_scale] <- tr$edge.length[edges_to_scale]*0.01

p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    mutate(highlight = ifelse(node %in% nodes_to_scale, T, F)) %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_color_manual(values = species_colors) +
    scale_x_continuous(limits = c(0, 0.0075)) +
    geom_treescale(x = .001, y = 15) +
    facet_grid2(~` `) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"),
        legend.position = "inside",
        legend.position.inside = c(.2, .8),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,-3,0,0), "mm")
    ) +
    guides(color = "none") +
    labs()

p1_1 <- isolates %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    select(genome_id, population, contig_species) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    ggplot() +
    geom_tile(aes(x = population, y = genome_id, fill = contig_species), color = "black", linewidth = .5) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors, breaks = c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.text = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.background = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
    labs()

# Panel B.


# Panel BC. Dxy ----
read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/gpa.csv"))
    gpar <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/gpar.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    sml <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/sml.csv"))
    spa <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/spa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/gene_order.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/gd.csv"))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)))
    gcn <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/gcn.csv"))
    cleaned_gene_names <- read_csv(paste0(folder_data, "genomics_analysis/species/gene_content/", set_name, "/cleaned_gene_names.csv"))
    return(list(gpa = gpa, gpar = gpar, list_sccg = list_sccg, sml = sml, spa = spa, gpatl = gpatl, gene_order = gene_order, gd = gd, gpacl = gpacl, gcn = gcn, cleaned_gene_names = cleaned_gene_names))
}
read_fsts <- function (set_name) {
    per_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/species/fst/", set_name,"/per_gene_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/species/fst/", set_name,"/per_locus_fst.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/species/fst/", set_name,"/gene_lengths.csv"))
    return(list(per_gene_fst = per_gene_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths))
}
read_dxys <- function (set_name) {
    gene_pop_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/dxy/", set_name,"/gene_pop_dxy.csv"))
    gene_ind_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/dxy/", set_name,"/gene_ind_dxy.csv"))
    return(list(gene_pop_dxy = gene_pop_dxy, gene_ind_dxy = gene_ind_dxy))
}
read_gcv_dxys <- function (set_name) {
    gcv_pop_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/gcv_dxy/", set_name, "/gcv_pop_dxy.csv"))
    gcv_ind_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/gcv_dxy/", set_name, "/gcv_ind_dxy.csv"))
    gcv_rep_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/gcv_dxy/", set_name, "/gcv_rep_dxy.csv"))
    return(list(gcv_pop_dxy = gcv_pop_dxy, gcv_ind_dxy = gcv_ind_dxy, gcv_rep_dxy = gcv_rep_dxy))
}
join_dists <- function (tt, dd, ff) {
    #' This function joins the genetic distance and geo distance
    gene_replicon <- filter(tt$gpacl, genome_id == tt$gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))

    dd$gene_ind_dxy %>%
        # dd$gene_ind_dxy %>%
        # left_join(dd$gene_lengths) %>%
        left_join(ff$gene_lengths) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
        left_join(sites_dist) %>%
        # left_join(sites_dist) %>%
        left_join(select(ff$per_gene_fst, gene, n_snps)) %>%
        mutate(dxy_scaled = dxy/sequence_length) %>%
        left_join(gene_replicon)
}
make_genome_wide_dxy <- function (xx) {
    #' This sums across the genes such that the final datables has choose(n,2) rows of genome-genome dxy
    xx %>%
        mutate(pops = paste0(population1, "-", population2)) %>%
        group_by(pops, genome_id1, genome_id2, dist_geo_km) %>%
        summarize(dxy = sum(dxy), sequence_length = sum(sequence_length), n_snps = sum(n_snps)) %>%
        mutate(dxy_scaled = dxy /sequence_length)
}
make_replicon_wide_dxy <- function (xx) {
    #' This sums across the genes such that the final datables has choose(n,2) rows of genome-genome dxy
    xx %>%
        mutate(pops = paste0(population1, "-", population2)) %>%
        group_by(pops, replicon_type, genome_id1, genome_id2, dist_geo_km) %>%
        summarize(dxy = sum(dxy), sequence_length = sum(sequence_length), n_snps = sum(n_snps)) %>%
        mutate(dxy_scaled = dxy /sequence_length)
}
make_dist_m <- function (xx_gn, dist_name) {
    #' Pivot the long format of distance to matrix format
    ml_gn <- xx_gn %>%
        ungroup() %>%
        select(genome_id1, genome_id2, {{dist_name}})

    m_gn <- ml_gn %>%
        bind_rows(rename(ml_gn, genome_id1 = genome_id2, genome_id2 = genome_id1)) %>%
        mutate(genome_id1 = factor(genome_id1, isolates$genome_id), genome_id2 = factor(genome_id2, isolates$genome_id)) %>%
        arrange(genome_id1, genome_id2) %>%
        pivot_wider(names_from = genome_id2, values_from = {{dist_name}}) %>%
        column_to_rownames(var = "genome_id1") %>%
        select(rownames(.)[1], everything()) %>%
        as.matrix

    diag(m_gn) <- 0
    return(m_gn)
}
plot_genome_wide_dxy <- function (xx_gn) {
    #' This plots the dxy
    xx_snps <- distinct(ungroup(xx_gn), n_snps)

    # Coerce the population groups
    xx_gn <- xx_gn %>%
        mutate(pops = case_when(
            pops == "low elevation-high elevation" ~ "high elevation-low elevation",
            pops == "urban-suburban" ~ "suburban-urban",
            T ~ pops
        )) %>%
        mutate(pops = str_remove_all(pops, " elevation"))

    # Mantel test
    m1 <- make_dist_m(xx_gn, dist_geo_km) # geo distance
    m2 <- make_dist_m(xx_gn, dxy_scaled) # genetic distance
    model <- mantel(m1, m2)
    r_squared <- model$statistic
    ast <- turn_p_to_asteriks(model$signif)

    xx_gn %>%
        ggplot() +
        geom_smooth(aes(x = dist_geo_km, y = dxy_scaled), method = "lm", se = F, color = "black") +
        #geom_point(aes(x = dist_geo_km, y = dxy_scaled, color = pops), shape = 21, size = 2, stroke = 1) +
        geom_point(aes(x = dist_geo_km, y = dxy_scaled), shape = 21, size = 2, stroke = .5) +
        #scale_color_manual(values = pops_colors) +
        #annotate("text", x = -Inf, y = Inf, label = paste(xx_snps$n_snps, "SNPs", ", r² =", round(r_squared, 3), ast), hjust = -.2, vjust = 1.5, size = 3, color = "black") +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            plot.title = element_text(size = 8)
        ) +
        guides() +
        labs(x = "Geographic distance (km)", y = "Dxy") +
        ggtitle(paste(xx_snps$n_snps, "SNPs", ", r² =", round(r_squared, 3), ast))
}
plot_replicon_wide_dxy <- function (xx_gn, xx_rp) {
    #' This plots the dxy per replicon
    xx_rp <- xx_rp %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        # Coerce the population groups
        mutate(pops = case_when(
            pops == "low elevation-high elevation" ~ "high elevation-low elevation",
            pops == "urban-suburban" ~ "suburban-urban",
            T ~ pops
        )) %>%
        mutate(pops = str_remove_all(pops, " elevation"))
    # Remove not assigned genes
    cat("\nDrop ", xx_rp$n_snps[which(is.na(xx_rp$replicon_type))[1]], " SNPs that do not belong any replicon")
    xx_rp <- drop_na(xx_rp, replicon_type)

    xx_snps <- distinct(ungroup(xx_rp), replicon_type, n_snps)

    # Mantel test
    m1 <- make_dist_m(xx_gn, dist_geo_km) # geo distance
    m2 <- make_dist_m(xx_gn, dxy_scaled) # genetic distance
    model <- mantel(m1, m2)
    r_squared <- model$statistic
    p_value <- turn_p_to_asteriks(model$signif)

    mantel_results <- xx_rp %>%
        group_by(replicon_type) %>%
        nest() %>%
        mutate(
            m1 = map(data, ~make_dist_m(.x, dist_geo_km)),
            m2 = map(data, ~make_dist_m(.x, dxy_scaled)),
            model = pmap(list(m1, m2), ~mantel(.x, .y)),
            r_squared = map_dbl(model, ~.x$statistic),
            ast = map_chr(model, ~turn_p_to_asteriks(.x$signif))
            # model = map(data, ~lm(dxy_scaled ~ dist_geo_km, data = .x)),
            # r_squared = map_dbl(model, ~ summary(.x)$r.squared),
            # p_value = map_dbl(model, ~ summary(.x)$coefficients[2, 4]),
            # ast = map_chr(p_value, ~ turn_p_to_asteriks(.x))
        ) %>%
        left_join(xx_snps)

    xx_rp %>%
        ggplot() +
        geom_smooth(aes(x = dist_geo_km, y = dxy_scaled), method = "lm", se = F, color = "black") +
        geom_point(aes(x = dist_geo_km, y = dxy_scaled), shape = 21, size = 2, stroke = .5) +
        geom_text(data = mantel_results, aes(label = paste(n_snps, "SNPs, r² =", round(r_squared, 3), ast)), x = -Inf, y = Inf, hjust = -.1, vjust = 1.5, size = 3, color = "black") +
        #scale_color_manual(values = pops_colors) +
        facet_grid2(~replicon_type) +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            legend.title = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            plot.title = element_text(size = 8),
            plot.background = element_blank()
        ) +
        guides() +
        labs(x = "Geographic distance (km)", y = "Dxy")
}
do_mantel <- function (xx_gn) {
    # Do mantel test
    m1 <- make_dist_m(xx_gn, dist_geo_km) # geo distance
    m2 <- make_dist_m(xx_gn, dxy_scaled) # genetic distance
    model <- mantel(m1, m2)
    return(tibble(replicon_type = "genome", r_squared = model$statistic, p_value = model$signif))
}
do_mantel_rep <- function (xx_rp) {
    # Do mantel test by replicon
    xx_rp %>%
        group_by(replicon_type) %>%
        nest() %>%
        mutate(
            m1 = map(data, ~make_dist_m(.x, dist_geo_km)),
            m2 = map(data, ~make_dist_m(.x, dxy_scaled)),
            model = pmap(list(m1, m2), ~mantel(.x, .y)),
            r_squared = map_dbl(model, ~.x$statistic),
            p_value = map_dbl(model, ~.x$signif)
        )
}

tb <- tibble(
    set_name = rep(c("elev_med", "urbn_mel"), each = 1)
) %>%
    mutate(
        tt = map(set_name, read_gpas),
        dd = map(set_name, read_dxys),
        ff = map(set_name, read_fsts),
        xx = pmap(list(tt, dd, ff), join_dists), # Per gene, dxy between two orthologs coming from two different genomes
        xx_gn = map(xx, make_genome_wide_dxy), # genome wide distance
        xx_rp = map(xx, ~make_replicon_wide_dxy(.x) %>% filter(replicon_type != "pAcce")),
        p_dxy = map(xx_gn, plot_genome_wide_dxy),
        p_dxy_rep = map2(xx_gn, xx_rp, plot_replicon_wide_dxy)
    )

# ----
p_right <- plot_grid(
    tb$p_dxy_rep[[2]], # meliloti
    tb$p_dxy_rep[[1]], # medicae
    nrow = 2, rel_widths = c(1,3), align = "h", axis = "tb",
    labels = LETTERS[2:5]
)

p <- plot_grid(
    plot_grid(p1, labels = "A"), p1_1 + guides(fill = "none"),
    p_right,
    #p2_1 + guides(fill = "none"), p2,
    nrow = 1,
    scale = .9, rel_widths = c(1,.1, 1.5)
    #align = "h", axis = "tb",
) +
    #draw_text("Single-copy core gene", x = .27, y = .95, size = 10, hjust = 0) +
    #draw_text("Gene content variation", x = .6, y = .95, size = 10, hjust = 0) +
    #draw_label("C", x = .79, y = .95, fontface = "bold") +
    draw_plot(get_legend(p1_1), x = -.4, y = .25) +
    #draw_plot(p3, width = .2, height = .35, x = .8, y = .6) +
    theme(plot.background = element_rect(color = NA, fill = "white"))


ggsave(here::here("plots/Fig4.png"), p, width = 12, height = 6)


#
tt <- read_gpas()
nrow(tt$gpa) # 26544

core <- tt$gpatl %>%
    group_by(gene) %>%
    filter(value == 1) %>%
    count() %>%
    ungroup() %>%
    filter(n == max(n))
nrow(core) / nrow(tt$gpa) *100

tt$gpacl %>%
    filter(str_detect(gene, "nod")) %>%
    filter(genome_id %in% paste0("g", c(2,3,15)))

