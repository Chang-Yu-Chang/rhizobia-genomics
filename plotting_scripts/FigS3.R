#' This script plots the pangenome composition: gene frequency spectrum and sampling regime
#' Within species

library(tidyverse)
library(cowplot)
library(ggsci)
source(here::here("metadata.R"))

iso <- read_csv(paste0(folder_data, "output/iso.csv"))
tt <- read_gpas()
gpatls <- tt$gpatl %>% left_join(select(iso, genome_id, contig_species))
compute_gfs <- function (gpatls) {
    gpatls %>%
        filter(value == 1) %>%
        group_by(gene) %>%
        count(name = "ngenomes") %>%
        group_by(ngenomes) %>%
        count() %>%
        mutate(n = n/1000)
}
count_genes <- function (gpatls) {
    gg <- gpatls %>%
        filter(value == 1) %>%
        group_by(gene) %>%
        count() %>%
        ungroup() %>%
        pull(n) %>%
        table()
    return(c(singleton = gg[[1]], core = last(gg), single_core = nrow(tt$list_sccg), total = sum(gg)))
}
plot_gfs <- function (gfs) {
    gfs %>%
        ggplot() +
        geom_hline(yintercept = 0) +
        geom_col(aes(x = ngenomes, y = n), fill = "white", color = "black", width = .8) +
        scale_x_continuous(breaks = seq(0, 40, 5)) +
        scale_y_continuous(breaks = c(0:10), limits = c(0, 6)) +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major.y = element_line(color = "grey90")
        ) +
        guides() +
        labs(x = "Number of genomes", y = "Number of genes (k)")

}


# Gene frequency spectrum  ----
# meliloti
p1_1 <- gpatls %>%
    filter(contig_species == "S. meliloti") %>%
    compute_gfs() %>%
    plot_gfs() +
    labs(title = "S. meliloti")

##
gpatls %>%
    filter(contig_species == "S. meliloti") %>%
    count_genes()

# medicae
p1_2 <- gpatls %>%
    filter(contig_species == "S. medicae") %>%
    compute_gfs() %>%
    plot_gfs() +
    labs(title = "S. medicae")

##
gpatls %>%
    filter(contig_species == "S. medicae") %>%
    count_genes()


# core vs accessory sampling ----
compute_pan <- function (mi) {
    # Single genome
    if(!is.matrix(mi)) return(tibble(core = sum(mi), total = sum(mi)))
    # Equal or more than two genomes
    ngenomes <- apply(mi, 2, sum)
    tng <- table(ngenomes)
    if (names(tng)[1] == "0") {
        tibble(core = last(tng)[1], total = sum(tng) - tng[1])
    } else if (names(tng)[1] == "1"){
        tibble(core = last(tng)[1], total = sum(tng))
    }
}
do_sampling <- function (gpa, n_boots = 100) {
    #tt <- read_gpas()
    m <- t(gpa[,-1])
    tb <- crossing(ngenome = 1:ncol(gpa[,-1]), boot = 1:n_boots) # Each sample of ngenome is repeated 100 times
    tb$pan <- NA

    for (i in 1:nrow(tb)) {
        set.seed(tb$boot[i])
        cat("", i)
        tb$pan[i] <- list(compute_pan(m[sample(1:nrow(m), tb$ngenome[i]),]))
    }

    tbp <- tb %>%
        unnest(cols = c(pan)) %>%
        group_by(ngenome) %>%
        reframe(enframe(quantile(core, c(0, 0.05, 0.5, 0.95, 1), na.rm = T), "quantile", "core"),
                enframe(quantile(total, c(0, 0.05, 0.5, 0.95, 1)), "quantile", "total")) %>%
        mutate(core = core/1000, total = total / 1000)

    tbpr <- tbp %>%
        filter(quantile %in% c("0%", "100%")) %>%
        pivot_longer(cols = c(core, total)) %>%
        pivot_wider(names_from = quantile, values_from = value)
    return(list(tbp = tbp, tbpr = tbpr))
}
plot_sampling <- function (tbp, tbpr) {
    tbp %>%
        group_by(ngenome) %>%
        pivot_longer(cols = c(-ngenome, -quantile)) %>%
        ggplot(aes(x = ngenome, y = value, color = name, linetype = quantile)) +
        geom_line() +
        geom_ribbon(data = tbpr, aes(x = ngenome, ymin = `0%`, ymax = `100%`, fill = name), inherit.aes = FALSE, alpha = 0.3) +
        scale_linetype_manual(values = c("0%" = 1, "5%" = 2, "50%" = 3, "95%" = 2, "100%" = 1)) +
        scale_x_continuous(breaks = seq(0, 40, 5)) +
        scale_y_continuous(breaks = seq(0, 30, 2), limits = c(4, 12)) +
        scale_color_aaas() +
        scale_fill_aaas() +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            legend.title = element_blank(),
            legend.position = "top",
            legend.margin = margin(0,0,0,0, "mm"),
            legend.box.margin = margin(0,0,0,0, "mm"),
            legend.background = element_rect(color = NA, fill = NA),
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_line(color = "grey90"),
            plot.margin = margin(0, 2,2,2, "mm")
        ) +
        guides(linetype = "none") +
        labs(x = "Number of genomes", y = "Number of genes (k)")
}

tb1 <- tt$gpa %>%
    select(gene, matches(paste0(iso$genome_id[iso$contig_species == "S. meliloti"], "$"))) %>%
    do_sampling(n_boots = 100)
p2_1 <- plot_sampling(tb1$tbp, tb1$tbpr)

tb2 <- tt$gpa %>%
    select(gene, matches(paste0(iso$genome_id[iso$contig_species == "S. medicae"], "$"))) %>%
    do_sampling(n_boots = 100)
p2_2 <- plot_sampling(tb2$tbp, tb2$tbpr)


p <- plot_grid(
    p1_1, p1_2, p2_1, p2_2 + theme(legend.position = "none"),
    nrow = 2, axis = "lrtb", align = "vh", scale = 0.9,
    labels = LETTERS[1:4], rel_widths = c(22,15)
) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS3.png"), p, width = 8, height = 6)


# Fraction of core ----
tb1$tbp$core[tb1$tbp$ngenome == max(tb1$tbp$ngenome) & tb1$tbp$quantile == "50%"] / tb1$tbp$total[tb1$tbp$ngenome == max(tb1$tbp$ngenome) & tb1$tbp$quantile == "50%"] * 100
tb2$tbp$core[tb2$tbp$ngenome == max(tb2$tbp$ngenome) & tb2$tbp$quantile == "50%"] / tb2$tbp$total[tb2$tbp$ngenome == max(tb2$tbp$ngenome) & tb2$tbp$quantile == "50%"] * 100
