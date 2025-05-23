#'

library(tidyverse)
library(ape)
library(TreeDist)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

# Prepare  data
tb <- tbtr %>%
    select(tr_type1 = tr_type, tr1 = tr) %>%
    cross_join(select(tbtr, tr_type2 = tr_type, tr2 = tr)) %>%
    mutate(tr_type1 = ordered(tr_type1, c("core", "gpa", "ani", "kmer")), tr_type2 = ordered(tr_type2, c("core", "gpa", "ani", "kmer"))) %>%
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

tb2 <- tb %>%
    mutate(rfd = map2_dbl(tr1, tr2, ~RobinsonFoulds(.x, .y, normalize = T, similarity = F))) %>%
    select(tr_type1, tr_type2, rfd, rfds) %>%
    unnest(rfds)

write_csv(tb2, paste0(folder_data, "phylogenomics_analysis/tree_distance/tb2.csv"))

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

p <- plot_hist_facet(tb2)
ggsave(paste0(folder_data, "phylogenomics_analysis/tree_distance/grfd_", set_name, ".png"), p, width = 5, height = 5)


