#' This script plots the figures of gene content

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

make_acce_fst <- function (per_acce_fst, gene_order, gpacl) {

    gene_replicon <- filter(gpacl, genome_id == gpacl$genome_id[1]) %>%
        select(gene, replicon_type) %>%
        replace_na(list(replicon_type = "others"))
    per_acce_fst %>%
        mutate(gene = factor(gene, gene_order$gene)) %>%
        left_join(gene_replicon) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce")))
}
plot_heatmap <- function (gpa, gpatl, gpacl, list_wgpa, by_replicon = F) {
    list_cg <- gpa$gene[apply(gpa[,-1], 1, sum) == ncol(gpa)-1]
    gene_order <- unique(gpacl$gene)
    p <- gpacl %>%
        #filter(!gene %in% list_cg) %>%
        # Exclude those genes that are assigned to different contigs in differnt genomes
        #filter(!gene %in% list_wgpa$gene) %>%
        left_join(isolates) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
        ggplot() +
        geom_tile(aes(x = gene, y = genome_id, fill = population)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_manual(values = population_colors) +
        facet_grid2(population~., scales = "free_y") +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 10, color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
        ) +
        guides(fill = "none") +
        labs(x = "gene cluster", y = "genome")

    if (by_replicon) {
        return(p + facet_grid(.~replicon_type, scales = "free_x", space = "free_x"))
    } else return(p)

}
plot_acce_fst <- function (acce_fst) {
    #' Plot accessory gene fst
    acce_fst %>%
        ggplot() +
        geom_segment(aes(x = gene, xend = gene, y = fst-0.01, yend = fst+0.01), alpha = 0.4) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
        ) +
        guides() +
        labs(x = "accessory gene", y = "Fst")
}
plot_gfs <- function (gpa) {
    #' Plot the gene frequency spectrum
    gfs <- gpa %>%
        pivot_longer(-gene) %>%
        group_by(gene) %>%
        summarize(n_genomes = sum(value)) %>%
        group_by(n_genomes) %>%
        count()
    gfs %>%
        ggplot() +
        geom_col(aes(x = n_genomes, y = n), color = "black", fill = "white") +
        scale_x_continuous(breaks = c(1, 31, seq(5, 30, 5))) +
        theme_classic() +
        theme() +
        guides() +
        labs(x = "# of genomes", y = "# of genes")
}
plot_singleton <- function (gpa, gpatl) {
    #' Plot the number of singletons per genome
    list_sg <- gpa$gene[which(apply(gpa[,-1], 1, sum) == 1)]
    gpatlsg <- gpatl %>%
        filter(gene %in% list_sg) %>%
        mutate(genome_id = factor(genome_id, isolates$genome_id))

    gpatlsg %>%
        group_by(genome_id) %>%
        count() %>%
        ggplot() +
        geom_col(aes(x = genome_id, y = n), color = "black") +
        theme_bw() +
        theme() +
        guides() +
        labs(title = "Symbiotic and non-symbiotic strains", y = "# of singltons")
}
plot_bcsm <- function (sml) {
    #' Plot the bray-curtis similarity matrix. Input is in long format
    sml %>%
        mutate(genome_id1 = factor(genome_id1, isolates$genome_id), genome_id2 = factor(genome_id2, rev(isolates$genome_id))) %>%
        ggplot() +
        geom_tile(aes(x = genome_id1, y = genome_id2, fill = bray_curtis_similarity)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradient(low = "snow", high = "maroon") +
        theme_bw() +
        theme() +
        labs()
}
get_gcn_agg <- function (gcn, cleaned_gene_names) {
    #' This aggregate the gene number for gene clusters that have the same gene symbol names
    gcn %>%
        left_join(cleaned_gene_names) %>%
        drop_na(from) %>%
        select(from, matches("g\\d")) %>%
        pivot_longer(-from) %>%
        mutate(name = factor(name, isolates$genome_id)) %>%
        group_by(from, name) %>%
        summarize(value = sum(value)) %>%
        pivot_wider() %>%
        ungroup()
}
plot_genecopy <- function (gcn_agg, labelled_genes) {
    #' This plots the heatmap of gene copy using the list of labelled genes
    tbw <- gcn_agg %>% filter(str_detect(from, labelled_genes))
    list_genes <- tbw$from
    tb <- tbw %>% pivot_longer(-from, names_to = "genome_id", values_to = "count")
    n_colors <- length(unique(tb$count))
    count_colors <- set_names(grey(seq(0.9, 0, length.out = n_colors)), sort(unique(tb$count)))

    tb %>%
        left_join(isolates) %>%
        # Orders
        mutate(from = factor(from, rev(list_genes))) %>%
        mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
        # Panels
        mutate(gene_family = str_sub(from, 1, 3)) %>%
        mutate(count = factor(as.character(count), 0:1000)) %>%
        ggplot() +
        geom_tile(aes(x = genome_id, y = from, fill = count), width = 0.95, height = 0.95) +
        scale_fill_manual(values = count_colors) +
        scale_x_discrete(position = "top", expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        facet_grid(gene_family~population, scale = "free", space = "free", switch = "y") +
        theme_minimal() +
        theme(
            #strip.text = element_text(angle = 90),
            strip.placement = "outside",
            strip.background = element_rect(color = "black"),
            panel.grid = element_blank(),
            panel.spacing = unit(3, "mm"),
            plot.background = element_rect(color = NA, fill = "white")
        ) +
        guides() +
        labs(x = "", y = "")
}
plot_bcsm_gcn <- function (gcn_agg) {
    #' Plot BC distance based on gene copy number
    nrow(gcn_agg) # 1990 genes
    m <- t(as.matrix(gcn_agg[,-1]))
    vegan::vegdist(m) %>%
        as.matrix() %>%
        as_tibble %>%
        mutate(genome_id1 = colnames(.)) %>%
        pivot_longer(-genome_id1, names_to = "genome_id2", values_to = "dist") %>%
        # left_join(rename_with(isolates, function (x) paste0(x, "1"))) %>%
        # left_join(rename_with(isolates, function (x) paste0(x, "2"))) %>%
        mutate(genome_id1 = factor(genome_id1, isolates$genome_id), genome_id2 = factor(genome_id2, rev(isolates$genome_id))) %>%
        ggplot() +
        geom_tile(aes(x = genome_id1, y = genome_id2, fill = dist)) +
        scale_x_discrete(position = "top",) +
        scale_fill_gradient(low = "white", high = "maroon", name = "BC distance") +
        theme_minimal() +
        theme(
            plot.background = element_rect(color = NA, fill = "white")
        ) +
        guides() +
        labs(x = "", y = "", title = paste0("copy number of ", nrow(gcn_agg), " genes"))
}
plot_wrapper <- function (set_name) {
    #set_name = "elev_med"
    #set_name = "urbn_mel"
    tt <- read_gpas(set_name)
    per_acce_fst <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/per_acce_fst.csv"))
    n_all <- nrow(tt$gene_order) # number of all genes
    n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes
    n_core = n_all-n_accessory
    acce_fst <- make_acce_fst(per_acce_fst, tt$gene_order, tt$gpacl)
    gcn_agg <- get_gcn_agg(tt$gcn, tt$cleaned_gene_names)

    p1 <- plot_heatmap(tt$gpa, tt$gpatl, tt$gpacl, list_wgpa) + ggtitle(paste0("Total: ", n_all, ", Core: ", n_core, ", Accessory: ", n_accessory))
    p2 <- plot_acce_fst(acce_fst)
    p <- plot_grid(p1, p2, nrow = 2, axis = "lrt", align = "vh")
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-01-gpa_heatmap.png"), p, width = 10, height = 6)
    p <- plot_gfs(tt$gpa)
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-02-gene_frequency_spectrum.png"), p, width = 5, height = 4)
    p <- plot_singleton(tt$gpa, tt$gpatl)
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-03-singletons.png"), p, width = 10, height = 8)
    p <- plot_bcsm(tt$sml)
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-04-gpa_sm.png"), p, width = 10, height = 8)
    p <- plot_genecopy(gcn_agg, "fdx|fix|hmp|nap|nod|nif|noe|nos|syr")
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-05-gcn_sym.png"), p, width = 8, height = 12)
    p <- plot_bcsm_gcn(gcn_agg)
    ggsave(paste0(folder_data, "genomics_analysis/gene_content/", set_name,"-06-gcn_sm.png"), p, width = 9, height = 8)
}

plot_wrapper("elev_med")
plot_wrapper("urbn_mel")
