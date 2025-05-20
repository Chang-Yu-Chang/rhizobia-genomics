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
    #distinct(name) %>%

gens %>%
    filter(gene == "arg") %>% view

target_info <- gens %>%
    left_join(select(contigs, seq_id = contig_id, replicon_type)) %>%
    filter(replicon_type == "pSymA") %>%
    filter(gene == "arg") %>%
    select(file_id, target_start = start, target_strand = strand)

# It's a mess to orginet the genome....
gens_new <- gens %>%
    left_join(target_info, by = "file_id") %>%
    filter(!is.na(target_start)) %>%
    group_by(file_id) %>%
    left_join(select(contigs, seq_id = contig_id, replicon_type)) %>%
    filter(replicon_type == "pSymA") %>%
    # For each genome...
    mutate(
        genome_length = max(end),
        # If target gene is on minus strand, invert all coordinates
        orientation_flipped = (target_strand == "-"),

        # Reorient, reassign positions
        new_start = if_else(
            orientation_flipped,
            genome_length - end,  # inverted positions
            start - target_start   # shift so target is at 0
        ),
        new_end = if_else(
            orientation_flipped,
            genome_length - start,
            end - target_start
        ),

        # Wrap around for circular genome
        new_start = if_else(new_start < 0, new_start + genome_length, new_start),
        new_end = if_else(new_end < 0, new_end + genome_length, new_end),

        # Set strand: all target genes should be on + strand
        new_strand = if_else(orientation_flipped, "+", "+")
    ) %>%
    ungroup() %>%
    # For genes in flipped genomes, invert the strand of other genes
    mutate(
        # For flipped genomes, invert strand of all genes
        final_strand = if_else(orientation_flipped, if_else(strand == "+", "-", "+"), strand)
    ) %>%
    select(-start, -end, -strand) %>%
    # Keep original strand
    select(file_id, seq_id, start = new_start, end = new_end, strand = final_strand, everything())

# Plot ----
p <- gggenomes(
    seqs = seqs,
    genes = gens_new %>% filter(str_detect(name, "nif|nod|fix"))
) +
    geom_bin_label() +
    geom_seq(aes(color = replicon_type), linewidth = 1) +
    geom_gene(aes(fill = name)) +
    scale_color_manual(values = c(chromosome = "black", pSymA = "darkred", pSymB = "darkblue", others = "grey80"), name = "Replicon")
p
ggsave(paste0(folder_data, "phylogenomics_analysis/comparative/01-genomes.png"), p, width = 10, height = 10)

