#'

library(tidyverse)
library(ape)
source(here::here("metadata.R"))
set.seed(42)

iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    arrange(contig_species)

# msa
list_genomes_symbiotic <- read_lines(paste0(folder_data, "genomics_analysis/count_snps/list_genomes_symbiotic.txt"))  # List of target genes
list_snps_symbiotic <- read_table(paste0(folder_data, "genomics_analysis/count_snps/list_snps_symbiotic.txt"))  # List of target genes
list_snps_symbiotic$GeneName <- str_remove(list_snps_symbiotic$GeneName, ".aln")
list_snps_symbiotic <- arrange(list_snps_symbiotic, GeneName)

cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=12)) %>% setNames(c(0:10, ">10"))

snps_dist <- list()
seq_len <- list()

for (gene_name in list_snps_symbiotic$GeneName) {
    # Read alignment
    aln <- read.FASTA(paste0(folder_genomics, "pangenome/aligned_gene_sequences/", gene_name, ".aln.fas"))
    names(aln) <- str_remove(names(aln), ";\\w+")
    aln <- aln[names(aln) %in% list_genomes_symbiotic] # subset those that are symbiotic
    seq_length <- length(aln[[1]])

    # Compute distance
    aln_dist <- dist.dna(aln, model = "N") %>%
        as.matrix %>% as_tibble %>%
        mutate(genome_id1 = names(aln)) %>%
        pivot_longer(-genome_id1, names_to = "genome_id2") %>%
        left_join(select(iso, genome_id1 = genome_id, species1 = contig_species), by = join_by(genome_id1)) %>%
        left_join(select(iso, genome_id2 = genome_id, species2 = contig_species), by = join_by(genome_id2)) %>%
        mutate(
            genome_id1 = factor(genome_id1, iso$genome_id),
            genome_id2 = factor(genome_id2, iso$genome_id)
        ) %>%
        filter(species1 == species2)

    snps_dist[[gene_name]] <- aln_dist
    seq_len[[gene_name]] <- seq_length

    # # plot
    # p <- aln_dist %>%
    #     mutate(
    #         genome_id1 = factor(genome_id1, rev(iso$genome_id)),
    #         value = ifelse(value > 10, ">10", value),
    #         value = factor(value)
    #     ) %>%
    #     ggplot() +
    #     geom_tile(aes(x = genome_id2, y = genome_id1, fill = value)) +
    #     scale_x_discrete(position = "top") +
    #     scale_fill_manual(values = cc, breaks = c(0:10, ">10")) +
    #     coord_cartesian(clip = "off") +
    #     theme_bw() +
    #     theme(
    #         axis.text.x = element_text(angle = 45, hjust = 0)
    #     ) +
    #     guides() +
    #     labs(title = paste0(gene_name, ":\t", seq_length, "bp"))
    #
    # ggsave(paste0(folder_data, "genomics_analysis/count_snps/distance_symbiotic/", gene_name, ".png"), p, width = 7, height = 6)

    cat("\n", gene_name)
}

#
snpd <- snps_dist %>%
    bind_rows(.id = "gene")
    #group_by(genome_id1, genome_id2, species1, species2) %>%
    #summarize(value = sum(value))
seql <- tibble(gene = names(seq_len), seql = unlist(seq_len))


# Filter for genes that might be of interest
snpd %>%
    left_join(seql) %>%
    filter(!str_detect(gene, "group")) %>%
    filter(value > 10, seql < 2000) %>%
    filter(!genome_id1 %in% c("g22", "g23"), !genome_id2 %in% c("g22", "g23")) %>%
    distinct(gene)

snpd %>%
    left_join(seql) %>%
    filter(!str_detect(gene, "group")) %>%
    filter(species1 == "S. medicae", species2 == "S. medicae") %>%
    filter(value > 10, seql < 2000)
    #filter(genome_id1 %in% c("g10", "g20"), genome_id2 %in% c("g20", "g10")) %>


snpd %>%
    left_join(seql) %>%
    filter(!str_detect(gene, "group")) %>%
    filter(species1 == "S. medicae", species2 == "S. medicae") %>%
    filter(genome_id1 == "g4", genome_id2 == "g17") %>%
    view
    #filter(value > 10, seql < 2000)


# snpd %>%
#     left_join(seql) %>%
#     filter(!str_detect(gene, "group")) %>%
#     filter(value > 10, seql < 2000) %>%
#     filter(!genome_id1 %in% c("g22", "g23"), !genome_id2 %in% c("g22", "g23")) %>%
#     distinct(gene) %>%
#     mutate(comm = paste0("mv distance_symbiotic/", gene, ".png candidates/", gene, ".png")) %>%
#     pull(comm) %>%
#     paste(collapse = ";")

iso %>%  filter(population == "VA") %>%
    view
