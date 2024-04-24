#' This script plots the pangenome composition: gene frequency spectrum and sampling regime

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("metadata.R"))

gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
gpat <- gpat %>% filter(!genome_id %in% c("g20", "g28"))
gpatl <- gpat %>%
    pivot_longer(-genome_id, names_to = "gene") %>%
    filter(value == 1)
m <- as.matrix(gpat[,-1])
colnames(m) <- NULL; dim(m)


# Plot GFS
p1 <- gpatl %>%
    group_by(gene) %>%
    count() %>%
    ggplot() +
    geom_histogram(aes(x = n), color = "black", fill = "white") +
    scale_x_continuous(breaks = c(1, seq(5, 30, 5))) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "Number of genomes", y = "Number of genes")

mt <- table(apply(m, 2, sum))
sum(mt[!names(mt) %in% c(0)]) # 26676
mt[names(mt) == 30] # 989

# Plot core vs accessory sampling

ngenomes <- apply(m, 2, sum)

compute_pan <- function (mi) {
    ngenomes <- apply(mi, 2, sum)
    tng <- table(ngenomes)
    tibble(core = last(tng)[1], accessory = sum(tng) - last(tng)[1] - tng[1])
}

compute_pan(m[1:5,])
n_boots <- 100
tb <- crossing(ngenome = 2:nrow(m), boot = 1:n_boots)
tb$pan <- NA

for (i in 1:nrow(tb)) {
    set.seed(tb$boot[i])
    cat("", i)
    tb$pan[i] <- list(compute_pan(m[sample(1:nrow(m), tb$ngenome[i]),]))
}

tbp <- tb %>%
    unnest(cols = c(pan)) %>%
    group_by(ngenome) %>%
    reframe(enframe(quantile(core, c(0.05, 0.5, 0.95)), "quantile", "core"),
            enframe(quantile(accessory, c(0.05, 0.5, 0.95)), "quantile", "accessory"))

tbpr <- tbp %>%
    filter(quantile %in% c("5%", "95%")) %>%
    pivot_longer(cols = c(core, accessory)) %>%
    pivot_wider(names_from = quantile, values_from = value)


p2 <- tbp %>%
    group_by(ngenome) %>%
    pivot_longer(cols = c(-ngenome, -quantile)) %>%
    ggplot(aes(x = ngenome, y = value, color = name, linetype = quantile)) +
    geom_line() +
    geom_ribbon(data = tbpr, aes(x = ngenome, ymin = `5%`, ymax = `95%`, fill = name), inherit.aes = FALSE, alpha = 0.3) +
    scale_linetype_manual(values = c("95%" = 1, "50%" = 2, "5%" = 1)) +
    scale_x_continuous(breaks = c(2, seq(5, 30, 5))) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "Number of genomes", y = "Number of genes")


p <- plot_grid(p1, p2, nrow = 2, labels = LETTERS[1:2], axis = "lr", align = "v", scale = 0.9) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS6.png"), p, width = 6, height = 6)











