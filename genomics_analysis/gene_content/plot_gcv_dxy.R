#' This script plots GCV dxy

library(tidyverse)
library(cowplot)
library(vegan) # for mantel test
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- select(isolates, -exp_id, -genome_name)
sites_dist <- read_csv(paste0(folder_data, "phenotypes/sites/sites_dist.csv"))

read_dxys <- function (set_name) {
    per_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_gene_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))

    gene_pop_dxy <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_pop_dxy.csv"))
    gene_ind_dxy <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_ind_dxy.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))

    gcv_pop_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gcv_pop_dxy.csv"))
    gcv_ind_dxy <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gcv_ind_dxy.csv"))

    return(list(per_gene_fst = per_gene_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths,
                gene_pop_dxy = gene_pop_dxy, gene_ind_dxy = gene_ind_dxy,
                gcv_pop_dxy = gcv_pop_dxy, gcv_ind_dxy = gcv_ind_dxy))
}
join_dists <- function (ind_dxy, gene_length, per_gene_fst, gpacl, sites_d, isolates) {
    #' This function joins the genetic distance and geo distance
    gene_replicon <- filter(gpacl, genome_id == gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))

    ind_dxy %>%
        # dd$gene_ind_dxy %>%
        # left_join(dd$gene_lengths) %>%
        left_join(gene_length) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
        left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
        left_join(sites_d) %>%
        #left_join(rename(gcv_ind_dxy, gcv_dxy = dxy)) %>%
        # left_join(sites_dist) %>%
        left_join(select(per_gene_fst, gene, n_snps)) %>%
        mutate(dxy_scaled = dxy/sequence_length) %>%
        left_join(gene_replicon)
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
        geom_point(aes(x = dist_geo_km, y = dxy_scaled, color = pops), shape = 21, size = 2, stroke = 1) +
        scale_color_manual(values = pops_colors) +
        annotate("text", x = -Inf, y = Inf, label = paste(xx_snps$n_snps, "SNPs", ", r² =", round(r_squared, 3), ast), hjust = -.2, vjust = 1.5, size = 3, color = "black") +
        theme_bw() +
        theme(
            legend.position = "right",
            legend.title = element_blank()
        ) +
        guides() +
        labs(x = "Geographic distance (km)", y = "Dxy")
}

#set_name = "elev_med"
#set_name = "urbn_mel"
tt <- read_gpas(set_name)
dd <- read_dxys(set_name)
n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes

dists <- dd$gcv_ind_dxy %>%
    rename(gcv_dxy = dxy) %>%
    left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
    left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
    left_join(sites_dist) %>%
    mutate(gcv_dxy_scaled = gcv_dxy / n_accessory)

m1 <- make_dist_m(dists, dist_geo_km) # geo distance
m2 <- make_dist_m(dists, gcv_dxy_scaled) # genetic distance
model <- mantel(m1, m2)
r_squared <- model$statistic
ast <- turn_p_to_asteriks(model$signif)

p <- dists %>%
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

ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-07-gcv_dxy.png"), p, width = 5, height = 4)

# # xx <- join_dists(dd$gene_ind_dxy, dd$gene_lengths, dd$per_gene_fst, tt$gpacl, sites_dist, isolates) # Per gene, dxy between two orthologs coming from two different genomes
# # xx %>%
# #     select(genome_id1, genome_id2)
# # xx_rp <- make_replicon_wide_dxy(xx)
# # xx_gn <- make_genome_wide_dxy(xx) # Genome wide, dxy between two genomes
# # nrow(xx_gn) # choose(10,2) or choose(17,2)
#
# p <- plot_genome_wide_dxy(xx_gn) + ggtitle(paste0(set_name, ": ", length(unique(c(xx_gn$genome_id1, xx_gn$genome_id2))), " genomes, ", nrow(xx_gn), " pairs"))
# ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-02-genome_dxy.png"), p, width = 5, height = 4)
# p <- plot_replicon_wide_dxy(xx_rp) + ggtitle(paste0(set_name, ": ", length(unique(c(xx_gn$genome_id1, xx_gn$genome_id2))), " genomes, ", nrow(xx_gn), " pairs"))
# ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-03-replicon_dxy.png"), p, width = 10, height = 4)
#
