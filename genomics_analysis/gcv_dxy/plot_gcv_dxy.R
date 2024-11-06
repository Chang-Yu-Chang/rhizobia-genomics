#' This script plots GCV dxy

library(tidyverse)
library(cowplot)
library(vegan) # for mantel test
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- select(isolates, -exp_id, -genome_name)
sites_dist <- read_csv(paste0(folder_data, "phenotypes/sites/sites_dist.csv"))

read_fsts <- function (set_name) {
    per_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_gene_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))
    return(list(per_gene_fst = per_gene_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths))
}
read_gcv_dxys <- function (set_name) {
    gcv_pop_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/gcv_pop_dxy.csv"))
    gcv_ind_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/gcv_ind_dxy.csv"))
    gcv_rep_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/gcv_rep_dxy.csv"))
    return(list(gcv_pop_dxy = gcv_pop_dxy, gcv_ind_dxy = gcv_ind_dxy, gcv_rep_dxy = gcv_rep_dxy))
}
join_dists <- function (dd) {
    gcv_rep_dxy_long <- dd$gcv_rep_dxy %>%
        pivot_wider(names_from = replicon_type, values_from = dxy, names_prefix = "dxy_")
    dd$gcv_ind_dxy %>%
        rename(gcv_dxy = dxy) %>%
        left_join(gcv_rep_dxy_long) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
        left_join(sites_dist) %>%
        mutate(gcv_dxy_scaled = gcv_dxy / n_accessory)
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
turn_p_to_asteriks <- function (p_value) {
    if (p_value < 0.001) {
        asterisks <- "***"
    } else if (p_value < 0.01) {
        asterisks <- "**"
    } else if (p_value < 0.05) {
        asterisks <- "*"
    } else {
        asterisks <- "ns"  # Not significant
    }
    return(asterisks)
}
plot_gcv_dxy <- function (dists) {
    dists %>%
        mutate(pops = paste0(population1, "-", population2)) %>%
        mutate(pops = case_when(
            pops == "low elevation-high elevation" ~ "high elevation-low elevation",
            pops == "urban-suburban" ~ "suburban-urban",
            T ~ pops
        )) %>%
        mutate(pops = str_remove_all(pops, " elevation")) %>%
        ggplot() +
        geom_smooth(aes(x = dist_geo_km, y = gcv_dxy_scaled), method = "lm", se = F, color = "black") +
        geom_point(aes(x = dist_geo_km, y = gcv_dxy_scaled, color = pops), shape = 21, size = 2, stroke = 1) +
        scale_color_manual(values = pops_colors) +
        annotate("text", x = -Inf, y = Inf, label = paste(n_accessory, " accessory genes", ", r² =", round(r_squared, 3), ast), hjust = -.2, vjust = 1.5, size = 3, color = "black") +
        theme_bw() +
        theme(
            legend.position = "right",
            legend.title = element_blank()
        ) +
        guides() +
        labs(x = "Geographic distance (km)", y = "GCV Dxy") +
        ggtitle(paste0(set_name, ": gene content, ", length(unique(c(dists$genome_id1, dists$genome_id2))), " genomes, ", nrow(dists), " pairs"))

}


#set_name = "elev_med"
#set_name = "urbn_mel"
tt <- read_gpas(set_name)
dd <- read_gcv_dxys(set_name)
n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes
n_accessory_rep <- tt$gpacl %>%
    select(replicon_type, gene) %>%
    distinct(replicon_type, gene) %>%
    group_by(replicon_type) %>%
    count()
dists <- join_dists(dd)

# Mantel test
m1 <- make_dist_m(dists, dist_geo_km) # geo distance
m2 <- make_dist_m(dists, gcv_dxy_scaled) # genetic distance
model <- mantel(m1, m2)
r_squared <- model$statistic
ast <- turn_p_to_asteriks(model$signif)

# Plot
p <- plot_gcv_dxy(dists)
ggsave(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name,"-01-gcv_dxy.png"), p, width = 5, height = 4)



dists %>%
    pivot_longer(cols = starts_with("dxy_"), names_prefix = "dxy_", names_to = "replicon_type", values_to = "dxy") %>%
    left_join(n_accessory_rep) %>%
    mutate(dxy_scaled = dxy/n) %>%
    select(replicon_type, genome_id1, genome_id2, dist_geo_km, dxy_scaled) %>%
    ggplot() +
    geom_smooth(aes(x = dist_geo_km, y = dxy_scaled), method = "lm") +
    geom_point(aes(x = dist_geo_km, y = dxy_scaled)) +
    facet_grid(~replicon_type) +
    theme_bw() +
    theme() +
    guides() +
    labs()


do_mantel <- function (data) {
    m1 <- make_dist_m(data, dist_geo_km) # geo distance
    m2 <- make_dist_m(data, dxy_scaled) # genetic distance
    model <- mantel(m1, m2)
    r_squared <- model$statistic
    ast <- turn_p_to_asteriks(model$signif)
    return(list(model, r_squared, ast))
}


xx <- dists %>%
    pivot_longer(cols = starts_with("dxy_"), names_prefix = "dxy_", names_to = "replicon_type", values_to = "dxy") %>%
    left_join(n_accessory_rep) %>%
    mutate(dxy_scaled = dxy/n) %>%
    select(replicon_type, genome_id1, genome_id2, dist_geo_km, dxy_scaled) %>%
    nest(data = -replicon_type) %>%
    #unnest(data) %>% view
    mutate(
        mod = map(data, do_mantel)
    )
xx$mod[[2]]



