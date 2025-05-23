#'

library(tidyverse)
library(cowplot)
library(ggh4x)
library(gggenomes)
source(here::here("metadata.R"))


load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
seqs <- read_csv(paste0(folder_data, "phylogenomics_analysis/comparative/seqs.csv"))
gens <- read_csv(paste0(folder_data, "phylogenomics_analysis/comparative/gens.csv"))
gens_oriented <- read_csv(paste0(folder_data, "phylogenomics_analysis/comparative/gens_oriented.csv"))

genome_ids <- c("g43", "g37", "g33", "g42", "g27", "g26", "g25", "g36", "g24", "g10", "g20", "g23", "g22", "g41", "g32", "g45", "g31", "g34", "g35", "g44", "g39", "g21", "g19", "g9", "g17", "g16", "g5", "g4", "g8", "g13", "g11", "g6", "g30", "g29", "g40", "g3", "g2", "g15")

# Panel A. contigs ----
p1 <- gggenomes(
    seqs = mutate(seqs, bin_id = factor(bin_id, genome_ids)) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "others"))),
    genes = gens_oriented %>% filter(str_detect(name, "nif|nod|fix"))
) +
    geom_bin_label() +
    geom_seq(aes(color = replicon_type), linewidth = 1) +
    geom_gene(aes(fill = name)) +
    #geom_gene_tag(aes(label = name), check_overlap = T) +
    scale_color_manual(values = c(chromosome = "black", pSymA = "darkred", pSymB = "darkblue", others = "grey80"), name = "Replicon") +
    #scale_x_continuous(breaks = paste0(seq(0,10, 2), "M")) +
    scale_y_continuous(expand = c(0,0), limits = c(.5, 38.5)) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.position = "top",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    ) +
    guides(fill = "none") +
    labs()
p1
# Panel B.  heatmap ----

iso <- read_csv(paste0(folder_data, "output/iso.csv"))
tt <- read_gpas()

genome_ids <- c("g43", "g37", "g33", "g42", "g27", "g26", "g25", "g36", "g24", "g10", "g20", "g23", "g22", "g41", "g32", "g45", "g31", "g34", "g35", "g44", "g39", "g21", "g19", "g9", "g17", "g16", "g5", "g4", "g8", "g13", "g11", "g6", "g30", "g29", "g40", "g3", "g2", "g15")

p2 <- tt$gd %>%
    filter(!str_detect(gene, "group")) %>%
    filter(str_detect(gene, "nod|nif|fix")) %>%
    mutate(ge = str_sub(gene, 1, 5) %>% str_remove("\\d$|_$")) %>%
    mutate(g = str_sub(ge, 1, 3)) %>%
    select(g, ge, gene, genome_id) %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    mutate(
        genome_id = factor(genome_id, rev(genome_ids)),
        contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens", "control"))
    ) %>%
    group_by(g, contig_species, genome_id, ge) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = ge, y = genome_id, fill = n)) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "grey80", high = "grey20") +
    facet_grid2(contig_species ~ g, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        legend.position = "right",
        strip.placement = "outside"
    ) +
    guides() +
    labs()


# ----
p <- plot_grid(
    p1, p2, nrow = 1,
    align = "h", axis = "tb",
    labels = LETTERS[1:2]
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig4.png"), p, width = 12, height = 6)

