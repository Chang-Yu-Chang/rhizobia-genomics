#' This script plots the pangenome composition: gene frequency spectrum and sampling regime

library(tidyverse)
library(cowplot)
library(ggsci)
source(here::here("metadata.R"))

iso <- read_csv(paste0(folder_data, "output/iso.csv"))

# Panel A. Gene frequency spectrum ----
tt <- read_gpas()

p1 <- tt$gpatl %>%
    filter(value == 1) %>%
    group_by(gene) %>%
    count(name = "ngenomes") %>%
    group_by(ngenomes) %>%
    count() %>%
    mutate(n = n/1000) %>%
    ggplot() +
    geom_hline(yintercept = 0) +
    geom_col(aes(x = ngenomes, y = n), fill = "white", color = "black", width = .8) +
    scale_x_continuous(breaks = seq(0, 40, 5)) +
    scale_y_continuous(breaks = c(0:10), limits = c(0, 9)) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey90")
    ) +
    guides() +
    labs(x = "Number of genomes", y = "Number of genes (k)")

# Number of core genes
gg <- table(apply(tt$gpa[,-1], 1, sum))
c(singleton = gg[[1]], core = last(gg), single_core = nrow(tt$list_sccg), total = sum(gg)) # 8875         987         775       26544

# Panel B. core vs accessory sampling ----
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
do_sampling <- function (n_boots = 100) {
    tt <- read_gpas()
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
        scale_x_continuous(breaks = seq(0, 40, 5)) +
        scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 27)) +
        scale_color_aaas() +
        scale_fill_aaas() +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            legend.title = element_blank(),
            legend.position = "right",
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

tb <- do_sampling()
p2 <- plot_sampling(tb$tbp, tb$tbpr)

p <- plot_grid(
    p1, p2,
    nrow = 2, axis = "lr", align = "v", scale = 0.9,
    labels = LETTERS[1:2]
) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS2.png"), p, width = 8, height = 6)
