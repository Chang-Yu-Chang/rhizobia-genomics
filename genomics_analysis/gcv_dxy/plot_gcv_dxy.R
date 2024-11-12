#' This script plots GCV dxy

library(tidyverse)
library(cowplot)
library(ggh4x)
library(vegan) # for mantel test
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- select(isolates, -exp_id, -genome_name)
sites_dist <- read_csv(paste0(folder_data, "phenotypes/sites/sites_dist.csv"))

read_gcv_dxys <- function (set_name) {
    pop_gcv_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/pop_gcv_dxy.csv"))
    ind_gcv_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/ind_gcv_dxy.csv"))
    rep_gcv_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/rep_gcv_dxy.csv"))
    n_acce <- read_csv(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name, "/n_acce.csv"))
    return(list(pop_gcv_dxy = pop_gcv_dxy, ind_gcv_dxy = ind_gcv_dxy, rep_gcv_dxy = rep_gcv_dxy, n_acce = n_acce))
}
join_gcv_dists <- function (dd) {
    # Replicons
    dd$rep_gcv_dxy %>%
        left_join(dd$n_acce) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
        #pivot_wider(names_from = replicon_type, values_from = gcv_dxy, names_prefix = "gcv_dxy_") %>%
        left_join(sites_dist) %>%
        mutate(gcv_dxy_scaled = gcv_dxy/n_accessory)

    # # Genome
    # dd$ind_gcv_dxy %>%
    #     rename(gcv_dxy_genome = gcv_dxy) %>%
    #     left_join(rep_gcv_dxy_long) %>%
    #     left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
    #     left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
    #     left_join(sites_dist) %>%
    #     mutate(gcv_dxy_scaled = gcv_dxy / n_accessory)
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
        geom_point(aes(x = dist_geo_km, y = gcv_dxy_scaled, color = pops), shape = 21, size = 2, stroke = 1) +
        scale_color_manual(values = pops_colors) +
        geom_text(data = mantel_results, aes(label = paste(n_accessory, "genes, r² =", round(r_squared, 3), ast)), x = -Inf, y = Inf, hjust = -.1, vjust = 1.5, size = 3, color = "black") +
        facet_grid2(~replicon_type) +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            plot.title = element_text(size = 8)
        ) +
        guides() +
        labs(x = "Geographic distance (km)", y = "GCV Dxy")
        #ggtitle(paste(n_accessory, " accessory genes", ", r² =", round(r_squared, 3), ast))
}

#set_name = "elev_med"
#set_name = "urbn_mel"
tt <- read_gpas(set_name)
dd <- read_gcv_dxys(set_name)
dists <- join_gcv_dists(dd)
#do_mantel(dists, gcv_dxy_scaled)

p <- plot_replicon_gcv_dxy(dists)
ggsave(paste0(folder_data, "genomics_analysis/gcv_dxy/", set_name,"-01-replicon_gcv_dxy.png"), p, width = 10, height = 3)


# dists %>%
#     pivot_longer(cols = starts_with("gcv_dxy_"), names_prefix = "gcv_dxy_", names_to = "replicon_type", values_to = "gcv_dxy") %>%
#     clean_pop_names() %>%
#     left_join(n_accessory_rep) %>%
#     mutate(dxy_scaled = dxy/n) %>%
#     #select(replicon_type, genome_id1, genome_id2, dist_geo_km, dxy_scaled) %>%
#     ggplot() +
#     geom_smooth(aes(x = dist_geo_km, y = dxy_scaled), method = "lm", se = F, color = "black") +
#     geom_point(aes(x = dist_geo_km, y = dxy_scaled, color = pops), shape = 21, size = 2, stroke = 1) +
#     scale_color_manual(values = pops_colors) +
#     facet_grid(~replicon_type) +
#     coord_cartesian(clip = "off") +
#     theme_classic() +
#     theme(
#         panel.border = element_rect(color = "black", fill = NA),
#         strip.background = element_blank()
#     ) +
#     guides() +
#     labs(x = "Geographic distance (km)", y = "GCV Dxy")
#
#
#
#
# xx <- dists %>%
#     pivot_longer(cols = starts_with("dxy_"), names_prefix = "dxy_", names_to = "replicon_type", values_to = "dxy") %>%
#     left_join(n_accessory_rep) %>%
#     mutate(dxy_scaled = dxy/n) %>%
#     select(replicon_type, genome_id1, genome_id2, dist_geo_km, dxy_scaled) %>%
#     nest(data = -replicon_type) %>%
#     #unnest(data) %>% view
#     mutate(
#         mod = map(data, do_mantel)
#     )
# xx$mod[[2]]
#
#
#
#
# # n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes
# # n_accessory_rep <- tt$gpacl %>%
# #     select(replicon_type, gene) %>%
# #     distinct(replicon_type, gene) %>%
# #     group_by(replicon_type) %>%
# #     count()
