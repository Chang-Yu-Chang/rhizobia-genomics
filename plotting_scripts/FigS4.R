#' This script plots the pangenome composition: gene frequency spectrum and sampling regime

library(tidyverse)
library(cowplot)
library(ggsci)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

# Gene frequency spectrum
plot_gfs <- function (set_name) {
    # Plot gene frequency spectrum
    tt <- read_gpas(set_name)
    tt$gpatl %>%
        filter(value == 1) %>%
        group_by(gene) %>%
        count(name = "ngenomes") %>%
        group_by(ngenomes) %>%
        count() %>%
        mutate(n = n/1000) %>%
        ggplot() +
        geom_hline(yintercept = 0) +
        geom_col(aes(x = ngenomes, y = n), fill = "white", color = "black", width = .8) +
        scale_x_continuous(breaks = c(1:ncol(tt$gpa[,-1]))) +
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

p_gfs1 <- plot_gfs("elev_med") + ggtitle("Elevation")
p_gfs2 <- plot_gfs("urbn_mel") + ggtitle("Urbanization")


count_genes <- function (set_name) {
    tt <- read_gpas(set_name)
    gg <- table(apply(tt$gpa[,-1], 1, sum))
    return(list(singleton = gg[1], core = last(gg), total = sum(gg)))
}

count_genes("elev_med")
count_genes("urbn_mel")

# Plot core vs accessory sampling
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
do_sampling <- function (set_name, n_boots = 100) {
    tt <- read_gpas(set_name)
    m <- t(tt$gpa[,-1])
    tb <- crossing(ngenome = 1:ncol(tt$gpa[,-1]), boot = 1:n_boots) # Each sample of ngenome is repeated 100 times
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
        scale_x_continuous(breaks = 1:20) +
        scale_y_continuous(breaks = seq(2, 12, 2), limits = c(5, 11)) +
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

tb1 <- do_sampling("elev_med")
tb2 <- do_sampling("urbn_mel")

p_s1 <- plot_sampling(tb1$tbp, tb1$tbpr)
p_s2 <- plot_sampling(tb2$tbp, tb2$tbpr)

p <- plot_grid(
    p_gfs1, p_gfs2, p_s1 + theme(legend.position = "none"), p_s2,
    nrow = 2, axis = "lrtb", align = "hv", scale = 0.9,
    labels = LETTERS[1:4], rel_widths = c(10,16)
) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS4.png"), p, width = 8, height = 6)
