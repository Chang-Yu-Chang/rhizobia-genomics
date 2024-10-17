#' This script plot the dxy per gene between 1) populations and 2) among genomes

library(tidyverse)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates <- select(isolates, -exp_id, -genome_name)
sites_dist <- read_csv(paste0(folder_data, "phenotypes/sites/sites_dist.csv"))
sites_dist <- mutate(sites_dist, dist_geo_km = dist_geo_m / 1000)

read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/spa.csv"))

    return(list(gpa = gpa, gene_order = gene_order, gpatl = gpatl, gpacl = gpacl, gd = gd, sml = sml, list_sccg = list_sccg, spa = spa))
}
read_dxys <- function (set_name) {
    gene_wide_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_wide_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))

    gene_pop_dxy <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_pop_dxy.csv"))
    gene_ind_dxy <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_ind_dxy.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))
    return(list(gene_wide_fst = gene_wide_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths,
                gene_pop_dxy = gene_pop_dxy, gene_ind_dxy = gene_ind_dxy))
}
join_dists <- function (ind_dxy, gene_length, gpacl, sites_d, isolates) {
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
        # left_join(sites_dist) %>%
        mutate(dxy_scaled = dxy/sequence_length) %>%
        left_join(gene_replicon)
}
make_replicon_wide_dxy <- function (xx) {
    #' This sums across the genes such that the final datables has choose(n,2) rows of genome-genome dxy
    xx %>%
        group_by(replicon_type, genome_id1, genome_id2, dist_geo_km) %>%
        summarize(dxy = sum(dxy), sequence_length = sum(sequence_length)) %>%
        mutate(dxy_scaled = dxy /sequence_length)
}
make_genome_wide_dxy <- function (xx) {
    #' This sums across the genes such that the final datables has choose(n,2) rows of genome-genome dxy
    xx %>%
        group_by(genome_id1, genome_id2, dist_geo_km) %>%
        summarize(dxy = sum(dxy), sequence_length = sum(sequence_length)) %>%
        mutate(dxy_scaled = dxy /sequence_length)
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
plot_ind_dxy <- function (xx) {
    xx %>%
        ggplot() +
        geom_point(aes(x = dist_geo_km, y = dxy_scaled), shape = 21) +
        theme_bw() +
        theme() +
        guides() +
        labs()
}
plot_genome_wide_dxy <- function (xx_gn) {
    #' This plots the dxy

    # Compute r^2
    model <- lm(dxy_scaled ~ dist_geo_km, data = xx_gn)
    r_squared <- summary(model)$r.squared
    p_value <- summary(model)$coefficients[2, 4] %>% turn_p_to_asteriks()

    xx_gn %>%
        ggplot() +
        geom_smooth(aes(x = dist_geo_km, y = dxy_scaled), method = "lm", se = F, color = "black") +
        geom_point(aes(x = dist_geo_km, y = dxy_scaled), shape = 21, size = 2) +
        annotate("text", x = -Inf, y = Inf, label = paste("R² =", round(r_squared, 3), p_value), hjust = -.3, vjust = 1.5, size = 3, color = "black") +
        theme_bw() +
        theme() +
        guides() +
        labs(x = "Geographic distance (km)", y = "Dxy")
}
plot_replicon_wide_dxy <- function (xx_rp) {
    #' This plots the dxy
    xx_rp <- xx_rp %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "psymA like", "psymB like", "others")))

    # Compute r^2
    lm_results <- xx_rp %>%
        group_by(replicon_type) %>%
        nest() %>%
        mutate(
            model = map(data, ~lm(dxy_scaled ~ dist_geo_km, data = .x)),
            r_squared = map_dbl(model, ~ summary(.x)$r.squared),
            p_value = map_dbl(model, ~ summary(.x)$coefficients[2, 4]),
            ast = map_chr(p_value, ~ turn_p_to_asteriks(.x))
        )

    xx_rp %>%
        ggplot() +
        geom_smooth(aes(x = dist_geo_km, y = dxy_scaled), method = "lm", se = F, color = "black") +
        geom_point(aes(x = dist_geo_km, y = dxy_scaled), shape = 21, size = 2) +
        geom_text(data = lm_results, aes(label = paste("R² =", round(r_squared, 3), ast)), x = -Inf, y = Inf, hjust = -.3, vjust = 1.5, size = 3, color = "black") +
        # annotate("text", x = -Inf, y = Inf, label = paste("R² =", round(r_squared, 3), p_value), hjust = -.3, vjust = 1.5, size = 3, color = "black") +
        facet_grid(~replicon_type) +
        theme_bw() +
        theme() +
        guides() +
        labs(x = "Geographic distance (km)", y = "Dxy")
}

#set_name = "elev_med"
set_name = "urbn_mel"
tt <- read_gpas(set_name)
dd <- read_dxys(set_name)

xx <- join_dists(dd$gene_ind_dxy, dd$gene_lengths, tt$gpacl, sites_dist, isolates) # Per gene, dxy between two orthologs coming from two different genomes
xx_rp <- make_replicon_wide_dxy(xx)
xx_gn <- make_genome_wide_dxy(xx) # Genome wide, dxy between two genomes
nrow(xx_gn) # choose(10,2) or choose(17,2)

p <- plot_genome_wide_dxy(xx_gn) + ggtitle(paste0(set_name, ": ", length(unique(c(xx_gn$genome_id1, xx_gn$genome_id2))), " genomes, ", nrow(xx_gn), " pairs"))
ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-02-genome_dxy.png"), p, width = 5, height = 4)
p <- plot_replicon_wide_dxy(xx_rp) + ggtitle(paste0(set_name, ": ", length(unique(c(xx_gn$genome_id1, xx_gn$genome_id2))), " genomes, ", nrow(xx_gn), " pairs"))
ggsave(paste0(folder_data, "genomics_analysis/fst/", set_name,"-03-replicon_dxy.png"), p, width = 10, height = 4)
