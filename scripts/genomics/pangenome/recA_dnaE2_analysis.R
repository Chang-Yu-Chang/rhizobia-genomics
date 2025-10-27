# Analyze recA and dnaE2 duplicates:
#   1. Summarize within-genome identity
#   2. Visualize phylogenetic trees

library(tidyverse)
library(ape)
library(ggtree)
library(cowplot)
source(here::here("metadata.R"))

out_dir <- paste0(folder_genomics, "pangenome/recA_dnaE2/")

# ---------- helper function ----------
summarize_identity <- function(file, gene) {
    message("Processing ", gene, "...")
    df <- read_tsv(
        file,
        col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen"),
        show_col_types = FALSE,
        progress = FALSE
    ) %>%
        filter(qseqid != sseqid, !is.na(pident)) %>%
        mutate(
            genome1 = str_extract(qseqid, "(?<=genomes/)(g[0-9]+)(?=\\.f)"),
            genome2 = str_extract(sseqid, "(?<=genomes/)(g[0-9]+)(?=\\.f)")
        ) %>%
        filter(
            !genome1 %in% c("g2", "g3", "g15"),
            !genome2 %in% c("g2", "g3", "g15"),
            genome1 == genome2        # keep only within-genome comparisons
        ) %>%
        group_by(genome = genome1) %>%
        summarise(
            n_pairs       = n(),
            mean_identity = mean(pident, na.rm = TRUE),
            max_identity  = max(pident, na.rm = TRUE)
        ) %>%
        arrange(desc(mean_identity))

    out_file <- file.path(out_dir, paste0(gene, "_within_identity_summary.tsv"))
    write_tsv(df, out_file)
    message("  → Wrote summary: ", out_file)
    print(df, n = 10)
    invisible(df)
}

# ---------- summarize BLAST results ----------
recA_id  <- summarize_identity(paste0(out_dir, "recA_self.tsv"),  "recA")
dnaE2_id <- summarize_identity(paste0(out_dir, "dnaE2_self.tsv"), "dnaE2")


# visualize phylogenetic trees
clean_tiplabel <- function(tree, drop_genomes = c("g2", "g3", "g15")) {
    genome_ids <- str_extract(tree$tip.label, "g[0-9]+")
    drop_idx <- which(is.na(genome_ids) | genome_ids %in% drop_genomes)
    if (length(drop_idx) > 0) {
        tree <- drop.tip(tree, drop_idx)
        genome_ids <- genome_ids[-drop_idx]
    }
    paralog_ids <- ave(seq_along(genome_ids), genome_ids, FUN = seq_along)
    tree$tip.label <- paste0(genome_ids, "_", paralog_ids)
    return(tree)
}
plot_with_taxalink <- function(tree, gene_name) {
    # Extract genome IDs and paralog counts
    tree_df <- tibble(label = tree$tip.label,
                      genome = str_extract(label, "g[0-9]+")) %>%
        add_count(genome, name = "n_paralogs")

    # Base plot
    p <- ggtree(tree, layout = "slanted") %<+% tree_df +
        geom_tiplab(size = 2) +
        ggtitle(paste0(gene_name)) +
        theme_tree() +
        scale_x_continuous(expand = c(0, 1))

    # Add taxalink for genomes with multiple paralogs
    multi_genomes <- tree_df %>%
        filter(n_paralogs > 1) %>%
        distinct(genome) %>%
        pull(genome)

    for (g in multi_genomes) {
        tips <- tree_df %>% filter(genome == g) %>% pull(label)
        if (length(tips) >= 2) {
            for (i in 1:(length(tips) - 1)) {
                p <- p + geom_taxalink(
                    taxa1 = tips[i],
                    taxa2 = tips[i + 1],
                    color = "grey40",
                    alpha = 0.6,
                    offset = 0.001,
                    outward = FALSE,
                    arrow = arrow(length = unit(0, "mm"))
                )
            }
        }
    }

    return(p)
}

recA_tree  <- read.tree(paste0(out_dir, "seqs/recA_tree.treefile"))
dnaE2_tree <- read.tree(paste0(out_dir, "seqs/dnaE2_tree.treefile"))
recA_tree  <- clean_tiplabel(recA_tree)
dnaE2_tree <- clean_tiplabel(dnaE2_tree)
p1 <- plot_with_taxalink(recA_tree,  "recA")
p2 <- plot_with_taxalink(dnaE2_tree, "dnaE2")
p <- plot_grid(p1, p2)

ggsave(paste0(folder_genomics, "pangenome/recA_dnaE2/trees.png"), p, width = 12, height = 8)
