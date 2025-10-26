#!/usr/bin/env Rscript
# =========================================================
# Analyze recA and dnaE2 duplicates:
#   1. Summarize within-genome identity
#   2. Visualize phylogenetic trees
# =========================================================

library(tidyverse)
library(ape)
library(ggtree)
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


# ---------- visualize phylogenetic trees ----------
message("\nDrawing phylogenetic trees...")

recA_tree  <- read.tree(paste0(out_dir, "seqs/recA_tree.treefile"))
dnaE2_tree <- read.tree(paste0(out_dir, "seqs/dnaE2_tree.treefile"))

ggtree(recA_tree) +
    geom_tiplab(size = 2) +
    ggtitle("recA Phylogeny")  +
    theme_tree2()

ggtree(dnaE2_tree) +
    geom_tiplab(size = 2) +
    ggtitle("dnaE2 Phylogeny") +
    theme_tree2()
