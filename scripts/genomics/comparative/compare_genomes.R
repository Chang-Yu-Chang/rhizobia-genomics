#'

library(tidyverse)
library(gggenomes)
source(here::here("metadata.R"))
set.seed(42)

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
genomes <- read_csv(paste0(folder_data, "genomics_analysis/genomes/genomes.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/genomes/contigs.csv"))

# Read files ----
list_seqs <- list()
list_genes <- list()
for (i in 1:nrow(iso)) {
#for (i in 3:4) {
    genome_id <- iso$genome_id[i]
    list_seqs[[genome_id]] <- read_seqs(paste0(folder_genomics, "fasta/genomes/", genome_id,".fasta"))
    list_genes[[genome_id]] <- read_feats(paste0(folder_genomics, "gff/", genome_id,".gff"))
    cat("\n", i)
}

# Bind files ----
genome_ids <- c("g43", "g37", "g33", "g42", "g27", "g26", "g25", "g36", "g24", "g10", "g20", "g23", "g22", "g41", "g32", "g45", "g31", "g34", "g35", "g44", "g39", "g21", "g19", "g9", "g17", "g16", "g5", "g4", "g8", "g13", "g11", "g6", "g30", "g29", "g40", "g3", "g2", "g15")

# Sequence and contigs
seqs <- bind_rows(list_seqs) %>%
    mutate(
        bin_id = file_id,
        bin_id = factor(bin_id, iso$genome_id),
        seq_id = paste0(file_id, "_", seq_id)
    ) %>%
    left_join(select(contigs, seq_id = contig_id, replicon, replicon_type)) %>%
    replace_na(list(replicon_type = "others")) %>%
    mutate(
        bin_id = factor(bin_id, genome_ids),
        replicon_type = ifelse(replicon_type == "pAcce", "others", replicon_type),
        replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "others"))
    ) %>%
    arrange(bin_id, replicon_type, desc(length)) %>%
    left_join(select(iso, file_id = genome_id, contig_species))

# genes
gens <- bind_rows(list_genes) %>%
    mutate(seq_id = paste0(file_id, "_", seq_id)) %>%
    mutate(name = str_remove(name, "_\\d+"))



write_csv(seqs, paste0(folder_data, "phylogenomics_analysis/comparative/seqs.csv"))
write_csv(gens, paste0(folder_data, "phylogenomics_analysis/comparative/gens.csv"))

# Orient the genes ----
seqs <- read_csv(paste0(folder_data, "phylogenomics_analysis/comparative/seqs.csv"))
gens <- read_csv(paste0(folder_data, "phylogenomics_analysis/comparative/gens.csv"))

# single copy core genes on psymA shared by all symbiotic genomes
gens %>%
    select(file_id, seq_id, gene) %>%
    left_join(select(contigs, seq_id = contig_id, replicon, replicon_type)) %>%
    drop_na(replicon_type) %>%
    distinct(replicon_type, gene, file_id) %>%
    filter(replicon_type == "pSymA") %>%
    group_by(gene) %>%
    count %>%
    filter(n == 35) %>%
    view

orient_replicon <- function (gens, repl, ge) {
    # Target gene to set at position 0
    target_info <- gens %>%
        left_join(select(contigs, seq_id = contig_id, replicon_type)) %>%
        filter(replicon_type == repl) %>%
        group_by(file_id) %>%
        mutate(contig_len = max(end)) %>%
        filter(gene == ge) %>%
        select(file_id, target_start = start, target_end = end, target_strand = strand, contig_len)

    list_orient <- list()

    for (i in 1:nrow(target_info)) {
        cat(i)
        gensi <- gens %>%
            left_join(select(contigs, seq_id = contig_id, replicon_type), by = join_by(seq_id)) %>%
            filter(file_id == target_info$file_id[i], replicon_type == repl) %>%
            select(-replicon_type)

        if (target_info$target_strand[i] == "+") {
            list_orient[[i]] <- gensi %>%
                # Orientate
                mutate(
                    new_start = start - target_info$target_start[i],
                    new_end = end - target_info$target_start[i],
                    new_start = ifelse(new_start < 0, new_start + target_info$contig_len[i], new_start),
                    new_end = ifelse(new_end < 0, new_end + target_info$contig_len[i], new_end),
                    new_strand = strand
                ) %>%
                arrange(new_start) %>%
                select(-start, -end, -strand) %>%
                select(file_id, seq_id, start = new_start, end = new_end, strand = new_strand, everything())

        } else if (target_info$target_strand[i] == "-") {
            list_orient[[i]] <- gensi %>%
                # Flip
                mutate(
                    new_start = target_info$contig_len[i] - end,
                    new_end = target_info$contig_len[i] - start,
                    new_strand = ifelse(strand == "+", "-", "+")
                ) %>%
                # Orientate
                mutate(
                    new_start = new_start - (target_info$contig_len[i] - target_info$target_end[i]), # the new target locus is also moved
                    new_end = new_end - (target_info$contig_len[i] - target_info$target_end[i]),
                    new_start = ifelse(new_start < 0, new_start + target_info$contig_len[i], new_start),
                    new_end = ifelse(new_end < 0, new_end + target_info$contig_len[i], new_end)
                ) %>%
                arrange(new_start) %>%
                select(-start, -end, -strand) %>%
                select(file_id, seq_id, start = new_start, end = new_end, strand = new_strand, everything())
        }
    }


    return(bind_rows(list_orient))
}
gens_pa <- orient_replicon(gens, "pSymA", "nifD_1")
gens_oriented <- gens %>%
    left_join(select(contigs, seq_id = contig_id, replicon_type), by = join_by(seq_id)) %>%
    # # All nonchr and non psyma contigs
    filter(replicon_type != "pSymA" | is.na(replicon_type)) %>%
    bind_rows(gens_pa) %>%
    select(file_id, seq_id, start, end, strand, type, name, gene)

write_csv(gens_oriented, paste0(folder_data, "phylogenomics_analysis/comparative/gens_oriented.csv"))

if (F) {


target_info <- gens %>%
    left_join(select(contigs, seq_id = contig_id, replicon_type)) %>%
    filter(replicon_type == "pSymA") %>%
    group_by(file_id) %>%
    mutate(contig_len = max(end)) %>%
    filter(gene == "nifD_1") %>%
    select(file_id, target_start = start, target_end = end, target_strand = strand, contig_len)

list_orient <- list()

for (i in 1:nrow(target_info)) {
    cat(i)
    gensi <- gens %>%
        left_join(select(contigs, seq_id = contig_id, replicon_type), by = join_by(seq_id)) %>%
        filter(file_id == target_info$file_id[i], replicon_type == "pSymA") %>%
        select(-replicon_type)

    if (target_info$target_strand[i] == "+") {
        list_orient[[i]] <- gensi %>%
            # Orientate
            mutate(
                new_start = start - target_info$target_start[i],
                new_end = end - target_info$target_start[i],
                new_start = ifelse(new_start < 0, new_start + target_info$contig_len[i], new_start),
                new_end = ifelse(new_end < 0, new_end + target_info$contig_len[i], new_end),
                new_strand = strand
            ) %>%
            arrange(new_start) %>%
            select(-start, -end, -strand) %>%
            select(file_id, seq_id, start = new_start, end = new_end, strand = new_strand, everything())

    } else if (target_info$target_strand[i] == "-") {
        list_orient[[i]] <- gensi %>%
            # Flip
            mutate(
                new_start = target_info$contig_len[i] - end,
                new_end = target_info$contig_len[i] - start,
                new_strand = ifelse(strand == "+", "-", "+")
            ) %>%
            # Orientate
            mutate(
                new_start = new_start - (target_info$contig_len[i] - target_info$target_end[i]), # the new target locus is also moved
                new_end = new_end - (target_info$contig_len[i] - target_info$target_end[i]),
                new_start = ifelse(new_start < 0, new_start + target_info$contig_len[i], new_start),
                new_end = ifelse(new_end < 0, new_end + target_info$contig_len[i], new_end)
            ) %>%
            arrange(new_start) %>%
            select(-start, -end, -strand) %>%
            select(file_id, seq_id, start = new_start, end = new_end, strand = new_strand, everything())
    }
}
}


# Plot ----
p <- gggenomes(
    seqs = seqs,
    genes = gens_oriented %>% filter(str_detect(name, "nod|nif|fix"))
) +
    geom_bin_label() +
    geom_seq(aes(color = replicon_type), linewidth = 1) +
    geom_gene(aes(fill = name)) +
    geom_gene_tag(aes(label = name), nudge_y=0.1, check_overlap = T, size = 2) +
    scale_color_manual(values = c(chromosome = "black", pSymA = "darkred", pSymB = "darkblue", others = "grey80"), name = "Replicon") +
    guides() +
    labs()

ggsave(paste0(folder_data, "phylogenomics_analysis/comparative/01-genomes.png"), p, width = 10, height = 10)

