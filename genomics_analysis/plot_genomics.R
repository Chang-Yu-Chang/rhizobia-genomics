#' This script plots the figures

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

plot_heatmap <- function (tt, by_replicon = F) {
    #set_name <- "elev_med"
    #set_name <- "urbn_mel"
    tt <- read_gpas(set_name)

    #list_cg <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) == ncol(tt$gpa)-1]
    gene_order <- levels(tt$gpacl$gene)
    n_all <- nrow(tt$gene_order) # number of all genes
    n_accessory <- tt$gpa$gene[apply(tt$gpa[,-1], 1, sum) != ncol(tt$gpa)-1] %>% length # number of accessory genes
    n_core = n_all-n_accessory

    # GPA
    p1 <- tt$gpacl %>%
        left_join(isolates) %>%
        mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        drop_na(replicon_type) %>%
        mutate(gene = factor(gene, rev(gene_order))) %>%
        ggplot() +
        geom_tile(aes(x = genome_id, y = gene, fill = population)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_fill_manual(values = population_colors) +
        facet_grid2(replicon_type~population, scales = "free", space = "free_y", switch = "y") +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            strip.placement = "outside",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
        ) +
        guides(fill = "none") +
        labs(x = "gene cluster", y = "genome")
        #ggtitle(paste0("Total: ", n_all, ", Core: ", n_core, ", Accessory: ", n_accessory))

    tb_ngenomes <- tt$gpacl %>%
        group_by(replicon_type, gene) %>%
        count(name = "n_genomes") %>%
        ungroup() %>%
        mutate(is_core = ifelse(n_genomes == max(n_genomes), "core", "accessory")) %>%
        arrange(gene) %>%
        ungroup()


    p2 <- tb_ngenomes %>%
        group_by(replicon_type, gene) %>%
        mutate(n_genomes_forplot = map(n_genomes, ~1:.x)) %>%
        unnest(n_genomes_forplot) %>%
        ungroup() %>%
        mutate(gene = factor(gene, rev(gene_order))) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        drop_na(replicon_type) %>%
        mutate(`shared genomes` = "shared genomes") %>%
        ggplot() +
        geom_tile(aes(x = n_genomes_forplot, y = gene), fill = "black") +
        scale_x_continuous(breaks = c(1, 5, 10, 15, 20), expand = c(0,0)) +
        facet_grid2(replicon_type~`shared genomes`, scales = "free", space = "free_y", strip = strip_vanilla(clip = "off")) +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(angle = -45, hjust = 1),
            strip.text.y = element_blank(),
            strip.clip = "off",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
            plot.background = element_blank()
        ) +
        guides() +
        labs(x = "gene cluster", y = "# of genomes")

    table(tb_ngenomes$is_core)
    p3 <- tb_ngenomes %>%
        mutate(gene = factor(gene, rev(gene_order))) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        drop_na(replicon_type) %>%
        mutate(dum = 1) %>%
        mutate(iscore = "is core") %>%
        ggplot() +
        geom_tile(aes(x = dum, y = gene, fill = is_core)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_fill_manual(values = c(core = "black", accessory = "white")) +
        facet_grid2(replicon_type~iscore, scales = "free", space = "free_y", strip = strip_vanilla(clip = "off")) +
        theme_classic() +
        coord_cartesian(clip = "off") +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(angle = -45, hjust = 1),
            strip.text.y = element_blank(),
            strip.clip = "off",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
            plot.background = element_blank()
        ) +
        guides(color = "none") +
        labs(x = "gene cluster", y = "")


    # if (by_replicon) p1 <- p1 + facet_grid(.~replicon_type, scales = "free_x", space = "free_x")

    p <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1, .1, .1), align = "h", axis = "tb") +
        theme(plot.background = element_rect(color = NA, fill = "white"))

    return(p)
}

set_name <- "elev_med"
set_name <- "urbn_mel"
tt <- read_gpas(set_name)
p <- plot_heatmap(tt)
ggsave(paste0(folder_data, "genomics_analysis/plot_genomics/", set_name,"-gcv.png"), p, width = 6, height = 8)
