#'

library(tidyverse)
library(ape)
library(TreeDist)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

# ACross all genomes ----
# Prepare data
tb <- tbtr %>%
    select(tr_type1 = tr_type, tr1 = tr) %>%
    cross_join(select(tbtr, tr_type2 = tr_type, tr2 = tr)) %>%
    mutate(tr_type1 = ordered(tr_type1, c("core", "gpa", "spa", "ani", "kmer")), tr_type2 = ordered(tr_type2, c("core", "gpa", "spa", "ani", "kmer"))) %>%
    mutate(rfds = NA) %>%
    filter(tr_type1 >= tr_type2)

# Randomize
set.seed(1)
for (j in 1:nrow(tb)) {
    rfds <- rep(NA, 1000)
    tr_q <- tb$tr2[[j]] # query tree

    for (i in 1:1000) {
        tr_i <- rtree(n = Ntip(tr_q))
        tr_i$tip.label <- tr_q$tip.label
        rfds[i] <- RobinsonFoulds(tb$tr1[[j]], tr_i, normalize = T)
    }
    cat("\n", j)
    tb$rfds[j] <- list(rfds)
}

tbb <- tb %>%
    mutate(rfd = map2_dbl(tr1, tr2, ~RobinsonFoulds(.x, .y, normalize = T, similarity = F))) %>%
    select(tr_type1, tr_type2, rfd, rfds) %>%
    unnest(rfds)

write_csv(tbb, paste0(folder_data, "phylogenomics_analysis/tree_distance/tbb.csv"))

plot_hist_facet <- function (tb2) {
    tb_obv <- tb2 %>%
        select(-rfds) %>%
        distinct(tr_type1, tr_type2, .keep_all = T)

    tb2 %>%
        ggplot() +
        geom_histogram(aes(x = rfds), binwidth = .05, color = "black", fill = NA) +
        geom_vline(data = tb_obv, aes(xintercept = rfd), color = "red", linetype = 2) +
        facet_grid2(tr_type1 ~ tr_type2, switch = "y", render_empty = F) +
        scale_x_continuous(breaks = seq(0,1,.2)) +
        scale_y_continuous(position = "right") +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        guides() +
        labs()

}

p <- plot_hist_facet(tbb)
ggsave(paste0(folder_data, "phylogenomics_analysis/tree_distance/grfd.png"), p, width = 5, height = 5)


# Within species ----
drop.tip(tb$tr1[[1]], "g2")
iso <- read_csv(paste0(folder_data, "output/iso.csv"))

tbtr1 <- tbtr %>%
    mutate(
        sp = "S. meliloti",
        tr = map2(tr, sp, ~drop.tip(.x, iso$genome_id[!iso$contig_species == .y]))
    )

tbtr2 <- tbtr %>%
    mutate(
        sp = "S. medicae",
        tr = map2(tr, sp, ~drop.tip(.x, iso$genome_id[!iso$contig_species == .y]))
    )

make_tb <- function (tbtr) {
    tbtr %>%
        select(tr_type1 = tr_type, tr1 = tr) %>%
        cross_join(select(tbtr, tr_type2 = tr_type, tr2 = tr)) %>%
        mutate(tr_type1 = ordered(tr_type1, c("core", "gpa", "spa", "ani", "kmer")), tr_type2 = ordered(tr_type2, c("core", "gpa", "spa", "ani", "kmer"))) %>%
        mutate(rfds = NA) %>%
        filter(tr_type1 >= tr_type2)
}

tb1 <- make_tb(tbtr1)
tb2 <- make_tb(tbtr2)


# Randomize
set.seed(1)
randomize_tree <- function (tb, n_random = 1000) {
    for (j in 1:nrow(tb)) {
        rfds <- rep(NA, n_random)
        tr_q <- tb$tr2[[j]] # query tree

        for (i in 1:n_random) {
            tr_i <- rtree(n = Ntip(tr_q))
            tr_i$tip.label <- tr_q$tip.label
            rfds[i] <- RobinsonFoulds(tb$tr1[[j]], tr_i, normalize = T)
        }
        cat("\n", j)
        tb$rfds[j] <- list(rfds)
    }

    tbb <- tb %>%
        mutate(rfd = map2_dbl(tr1, tr2, ~RobinsonFoulds(.x, .y, normalize = T, similarity = F))) %>%
        select(tr_type1, tr_type2, rfd, rfds) %>%
        unnest(rfds)
    return(tbb)
}
tbb1 <- randomize_tree(tb1)
tbb2 <- randomize_tree(tb2)

p <- plot_hist_facet(tbb1) + labs(title = "S. meliloti")
ggsave(paste0(folder_data, "phylogenomics_analysis/tree_distance/grfd_meliloti.png"), p, width = 5, height = 5)
p <- plot_hist_facet(tbb2) + labs(title = "S. medicae")
ggsave(paste0(folder_data, "phylogenomics_analysis/tree_distance/grfd_medicae.png"), p, width = 5, height = 5)

