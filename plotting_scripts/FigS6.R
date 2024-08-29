#' This script plots the pangenome composition: gene frequency spectrum and sampling regime

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
library(ggsci)
source(here::here("metadata.R"))

gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpatl.csv"))
m <- as.matrix(gpat[,-1])
colnames(m) <- NULL; dim(m)

# Plot GFS
p1 <- gpatl %>%
    group_by(gene) %>%
    count(name = "ngenomes") %>%
    group_by(ngenomes) %>%
    count() %>%
    mutate(n = n/1000) %>%
    ggplot() +
    geom_col(aes(x = ngenomes, y = n), fill = "white", color = "black", width = 1) +
    scale_x_continuous(breaks = c(1, seq(5, 36, 5))) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs(x = "Number of genomes", y = "Number of genes (k)")

mt <- table(apply(m, 2, sum))
sum(mt[!names(mt) %in% c(0)]) # 26504 all genes
mt[names(mt) == length(unique(gpatl$genome_id))] # 988 core gene

# Plot core vs accessory sampling
ngenomes <- apply(m, 2, sum)
compute_pan <- function (mi) {
    # Single genome
    if(!is.matrix(mi)) return(tibble(core = NA, total = sum(mi)))
    # Equal or more than two genomes
    ngenomes <- apply(mi, 2, sum)
    tng <- table(ngenomes)
    if (names(tng)[1] == "0") {
        tibble(core = last(tng)[1], total = sum(tng) - tng[1])
    } else if (names(tng)[1] == "1"){
        tibble(core = last(tng)[1], total = sum(tng))
    }
}

#compute_pan(m[1,])
n_boots <- 100
tb <- crossing(ngenome = 1:nrow(m), boot = 1:n_boots) # Each sample of ngenome is repeated 100 times
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


p2 <- tbp %>%
    group_by(ngenome) %>%
    pivot_longer(cols = c(-ngenome, -quantile)) %>%
    ggplot(aes(x = ngenome, y = value, color = name, linetype = quantile)) +
    geom_line() +
    geom_ribbon(data = tbpr, aes(x = ngenome, ymin = `0%`, ymax = `100%`, fill = name), inherit.aes = FALSE, alpha = 0.3) +
    scale_linetype_manual(values = c("0%" = 1, "5%" = 2, "50%" = 3, "95%" = 2, "100%" = 1)) +
    scale_x_continuous(breaks = c(1, seq(5, 36, 5))) +
    scale_y_continuous(breaks = seq(0, 25, 5)) +
    scale_color_aaas() +
    scale_fill_aaas() +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.3),
        legend.background = element_rect(color = "black", fill = "white")
    ) +
    guides(linetype = "none") +
    labs(x = "Number of genomes", y = "Number of genes (k)")

p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], axis = "lr", align = "v", scale = 0.9) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS6.png"), p, width = 8, height = 4)
