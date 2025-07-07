#' Isolation by distance

library(tidyverse)
library(cowplot)
library(ggh4x)
library(vegan) # for mantel test
library(lme4) # for lmer
library(car) # for anova
library(broom.mixed)
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- select(isolates, -exp_id, -genome_name)
sites_dist <- read_csv(paste0(folder_data, "phenotypes/sites/sites_dist.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))

# Figures ----
# core ----
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
    pop_gcv_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/gcv_dxy/", set_name, "/pop_gcv_dxy.csv"))
    ind_gcv_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/gcv_dxy/", set_name, "/ind_gcv_dxy.csv"))
    rep_gcv_dxy <- read_csv(paste0(folder_data, "genomics_analysis/species/gcv_dxy/", set_name, "/rep_gcv_dxy.csv"))
    n_acce <- read_csv(paste0(folder_data, "genomics_analysis/species/gcv_dxy/", set_name, "/n_acce.csv"))
    return(list(pop_gcv_dxy = pop_gcv_dxy, ind_gcv_dxy = ind_gcv_dxy, rep_gcv_dxy = rep_gcv_dxy, n_acce = n_acce))
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

# GCV ----
join_gcv_dists <- function (dd) {
    # Replicons
    dd$rep_gcv_dxy %>%
        left_join(dd$n_acce) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
        #pivot_wider(names_from = replicon_type, values_from = gcv_dxy, names_prefix = "gcv_dxy_") %>%
        left_join(sites_dist) %>%
        mutate(gcv_dxy_scaled = gcv_dxy/n_accessory) %>%
        mutate(pops = paste0(population1, "-", population2)) %>%
        mutate(pops = case_when(
            pops == "low elevation-high elevation" ~ "high elevation-low elevation",
            pops == "urban-suburban" ~ "suburban-urban",
            T ~ pops
        )) %>%
        mutate(pops = str_remove_all(pops, " elevation"))
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
make_genome_wide_gcv_dxy <- function (dists) {
    dists %>%
        group_by(pops, genome_id1, genome_id2) %>%
        summarize(gcv_dxy_scaled = sum(gcv_dxy_scaled))
}
do_mantel <- function (data, genetic_d) {
    m1 <- make_dist_m(data, dist_geo_km) # geo distance
    m2 <- make_dist_m(data, {{genetic_d}}) # genetic distance
    model <- mantel(m1, m2)
    return(tibble(r_squared = model$statistic, p_value = model$signif))
}
clean_pop_names <- function (tb) {
    tb %>%
        mutate(pops = paste0(population1, "-", population2)) %>%
        mutate(pops = case_when(
            pops == "low elevation-high elevation" ~ "high elevation-low elevation",
            pops == "urban-suburban" ~ "suburban-urban",
            T ~ pops
        )) %>%
        mutate(pops = str_remove_all(pops, " elevation"))
}
plot_replicon_gcv_dxy <- function (dists) {
    mantel_results <- dists %>%
        group_by(replicon_type) %>%
        nest() %>%
        mutate(
            mod = map(data, ~do_mantel(.x, gcv_dxy_scaled))
            # r_squared = map_dbl(mod, ~.x$statistic),
            # ast = map_chr(mod, ~turn_p_to_asteriks(.x$signif))
            # model = map(data, ~lm(dxy_scaled ~ dist_geo_km, data = .x)),
            # r_squared = map_dbl(model, ~ summary(.x)$r.squared),
            # p_value = map_dbl(model, ~ summary(.x)$coefficients[2, 4]),
            # ast = map_chr(p_value, ~ turn_p_to_asteriks(.x))
        ) %>%
        unnest(mod) %>%
        mutate(ast = map_chr(p_value, turn_p_to_asteriks)) %>%
        left_join(distinct(dists, replicon_type, n_accessory)) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce")))


    dists %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        clean_pop_names() %>%
        ggplot() +
        geom_smooth(aes(x = dist_geo_km, y = gcv_dxy_scaled), method = "lm", se = F, color = "black") +
        geom_point(aes(x = dist_geo_km, y = gcv_dxy_scaled), shape = 21, size = 2, stroke = .5) +
        #scale_color_manual(values = pops_colors) +
        geom_text(data = mantel_results, aes(label = paste(n_accessory, "genes, r² =", round(r_squared, 3), ast)), x = -Inf, y = Inf, hjust = -.1, vjust = 1.5, size = 3, color = "black") +
        facet_grid2(~replicon_type) +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            plot.title = element_text(size = 8),
            plot.background = element_blank()
        ) +
        guides() +
        labs(x = "Geographic distance (km)", y = "GCV Dxy")
    #ggtitle(paste(n_accessory, " accessory genes", ", r² =", round(r_squared, 3), ast))
}

tbg <- tibble(set_name = c("elev_med", "urbn_mel")) %>%
    mutate(
        tt = map(set_name, read_gpas),
        dd = map(set_name, read_gcv_dxys),
        dists = map(dd, ~join_gcv_dists(.x) %>% filter(replicon_type != "pAcce")),
        xx_gn = map(dists, ~make_genome_wide_gcv_dxy(.x)),
        p = map(dists, plot_replicon_gcv_dxy)
    )


# ----

theme_consist <- function () {
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.02, .65),
        legend.justification.inside = "left",
        legend.text = element_text(size = 8),
        legend.background = element_rect(color = "black", fill = "snow", linewidth = .3),
        # legend.box.margin = margin(0,1,0,0, "mm"),
        # legend.margin = margin(0,1,0,0, "mm"),
        legend.key = element_blank(),
        plot.title = element_text(size = 15, face = "italic")
    )
}
p_combined <- plot_grid(
    tb$p_dxy_rep[[2]] + theme_consist() + xlim(0, 17) + ylim(0, 0.02) + ggtitle("S. meliloti"),
    tb$p_dxy_rep[[1]] + theme_consist() + xlim(0, 17) + ylim(0, 0.02) + ggtitle("S. medicae"),
    tbg$p[[2]] + theme_consist() + theme(legend.position = "none") + xlim(0, 17) + ylim(0, .7) + ggtitle("S. meliloti"),
    tbg$p[[1]] + xlim(0, 17) + ylim(0, .7) + theme(legend.position = "none", plot.title = element_text(size = 15, face = "italic")) + ggtitle("S. medicae"),
    scale = 0.95, ncol = 1, align = "hv", axis = "trl",
    rel_heights = c(1,1,1,1.2), labels = LETTERS[1:4], label_x = .05
)

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig4.png"), scale = 1) +
    draw_plot(p_combined, x = .05, width = .95) +
    theme(plot.background = element_rect(fill = "white", color = NA))


ggsave(here::here("plots/Fig4.png"), p, width = 8, height = 10)


# Stat ----
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

