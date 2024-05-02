#' This script plots the genomes

renv::load()
library(tidyverse)
library(gggenomes)
library(cowplot)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/synteny/gggenomes.rdata"))
# mss <- read_csv(paste0(folder_data, "phylogenomics_analysis/synteny/mss.csv"))
# mgs <- read_csv(paste0(folder_data, "phylogenomics_analysis/synteny/mgs.csv"))

unique(mgs$replicon)
ms <- mss %>%
    #filter(replicon == "psyma")
    #filter(str_detect(seq_id, "g5"))
    {.}
mg <- mgs %>%
    #filter(replicon == "psyma") %>%
    #filter(str_detect(seq_id, "g5")) %>%
    filter(str_detect(name, "nif|nod|fix")) %>%
    mutate(gene_group = str_sub(name, 1, 3))

p1 <- gggenomes(
    genes = filter(mg, replicon == "chromosome") %>% filter(!str_detect(seq_id, "g2_|g3_|g15_")),
    seqs = filter(ms, replicon == "chromosome") %>% filter(!str_detect(seq_id, "g2_|g3_|g15_"))
) +
    geom_seq() +
    geom_bin_label() +
    geom_gene(aes(color = gene_group)) +
    theme(legend.position = "none", panel.border = element_rect(color = 1, fill = NA)) +
    labs(title = "chromosome")
p2 <- gggenomes(
    genes = filter(mg, replicon == "psyma"),
    seqs = filter(ms, replicon == "psyma")
) +
    geom_seq() +
    geom_bin_label() +
    geom_gene(aes(color = gene_group)) +
    theme(legend.position = "bottom", panel.border = element_rect(color = 1, fill = NA)) +
    labs(title = "pSymA or alike")
p3 <- gggenomes(
    genes = filter(mg, replicon == "psymb"),
    seqs = filter(ms, replicon == "psymb")
) +
    geom_seq() +
    geom_bin_label() +
    geom_gene(aes(color = gene_group)) +
    theme(legend.position = "none", panel.border = element_rect(color = 1, fill = NA)) +
    labs(title = "pSymB or alike")

p <- plot_grid(p1, p2, p3, nrow = 1, axis = "tb", align = "h")

ggsave(paste0(folder_data, "phylogenomics_analysis/synteny/01-sym_genes.png"), p, width = 15, height = 10)

# 2. plot the sym genes in one ggplot


# mss <- bind_rows(tb$seqs_chr[1:3])
# mgs <- bind_rows(tb$genes_chr[1:3]) %>%
#     mutate(gene_group = str_sub(name, 1, 3)) %>%
#     #drop_na(gene_group)
#     filter(str_detect(gene_group, "rps|dna|rec|rpo"))

p <- gggenomes(genes = mg, seqs = ms) +
    geom_seq() +
    geom_bin_label() +
    geom_gene(aes(color = gene_group)) +
    facet_grid(.~replicon, scales = "free", space = "free_x") +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(title = "chromosome")
ggsave(paste0(folder_data, "phylogenomics_analysis/synteny/02-sym_genes.png"), p, width = 10, height = 10)

# 3. plot the PAP tree with pSymA
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
library(ggtree)
library(tidytree)
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
isolates_tax <- mutate(isolates_tax, contig_species = str_remove(contig_species, "E. "))
p1 <- tr_acce %>%
    drop.tip(c("g2", "g3", "g15")) %>%
    as_tibble() %>%
    left_join(rename(isolates_tax, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(hjust = 0, offset = 0.001) +
    geom_tippoint(aes(color = contig_species)) +
    # geom_label2(aes(subset=(node %in% c(32,34,35))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    # geom_label2(aes(subset=(node %in% c(32,34,35))), label = c("Ensifer meliloti", "Ensifer medicae", "Ensifer spp."), label.size = 0, fill = NA, nudge_x = -.02, nudge_y = 1, hjust = 1) +
    annotate("segment", x = 0.1, xend = 0.2, y = 3, yend = 3, arrow = arrow(angle = 90, ends = "both", length = unit(1, "mm"))) +
    annotate("text", x = 0.15, y = 2, label = "0.1") +
    scale_x_continuous(limits = c(0, 1.03), expand = c(0,.001)) +
    scale_y_continuous(limits = c(1, 27), expand = c(0,1)) +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = c(0.3, 0.9),
        legend.background = element_rect(color = "black", fill = "white"),
        axis.line.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.length = unit(0, "mm"),
        axis.text.y = element_text(size = 10)
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(title = "gene content tree")

ordered_tips <- p1$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)

p2 <- gggenomes(
    genes = filter(mg, replicon == "psyma"),
    seqs = filter(ms, replicon == "psyma") %>%
        filter(file_id != "g20") %>%
        mutate(file_id = factor(file_id, rev(ordered_tips))) %>% arrange(file_id)
) +
    geom_seq() +
    geom_bin_label() +
    geom_gene(aes(color = gene_group, fill = gene_group)) +
    scale_y_continuous(limits = c(1, 27), expand = c(0,1)) +
    theme(legend.position = "bottom") +
    labs(title = "pSymA or alike")

p <- plot_grid(p1, p2, nrow = 1, axis = "tb", align = "h")
ggsave(paste0(folder_data, "phylogenomics_analysis/synteny/03-tree_psyma.png"), p, width = 15, height = 10)















if (F) {

    # Bind gff into one table
    list_seqs <- rep(list(NA), nrow(isolates))
    list_genes <- rep(list(NA), nrow(isolates))

    for (i in c(1:22, 24:32)) {
        print(i)
        # Fasta
        m_seqs <- read_seq_len(paste0(folder_genomics, "fasta/genomes/", isolates$genome_id[i],".fasta")) %>%arrange(desc(length))
        #m_seqs_sub <- m_seqs %>% filter(length > 1000000)

        # GFF
        m_genes <- read_gff3(paste0(folder_genomics, "gff/genomes/", isolates$genome_id[i],".gff"))
        # m_genes_sub <- m_genes %>%
        #     filter(str_detect(gene, "nod|fix|nif")) %>%
        #     mutate(gene_type = str_sub(gene, 1, 3)) %>%
        #     {.}

        #
        list_seqs[[i]] <- mutate(m_seqs, genome_id = isolates$genome_id[i]) %>% select(genome_id, everything())
        list_genes[[i]] <- mutate(m_genes, genome_id = isolates$genome_id[i]) %>% select(genome_id, everything())

        # list_p[[i]] <- gggenomes(seqs = m_seqs_sub, genes = m_genes_sub) +
        #     geom_seq() + geom_bin_label() +
        #     geom_gene(aes(color = gene_type, fill = gene_type)) +
        #     theme()
    }

    bind_rows(list_seqs[!is.na(list_seqs)])

}

if (FALSE) {
    data(package="gggenomes")

    emale_seqs
    # emale_seqs <- emale_seqs %>%
    #   extract(seq_desc, into = c("emale_type", "is_typespecies"), "=(\\S+) \\S+=(\\S+)", remove=F, convert=T) %>%
    #   dplyr::arrange(emale_type, length)

    emale_genes <- emale_genes %>%
        mutate(gc_content = as.numeric(gc_content))
    p <- gggenomes(genes = emale_genes, seqs = emale_seqs, feats = emale_tirs, links = emale_ava) +
        geom_seq() +
        geom_bin_label() +
        geom_feat(size = 5) +
        geom_gene(aes(fill = gc_content)) +
        geom_link(color = NA) +
        scale_fill_distiller(palette = "Spectral") +
        theme()
    p
    p %>%
        add_feats(emale_gc) +
        geom_ribbon(aes(x=(x+xend)/2, ymax=y+.24, ymin=y+.38-(.4*score),group=seq_id,
                        linetype="GC-content"), feats(emale_gc),fill="blue", alpha=.5)
    emale_cogs
    emale_ava

    gggenomes(genes = emale_genes, seqs = emale_seqs, feats = emale_tirs, links = emale_ava) %>%
        add_feats(emale_gc) %>%
        add_clusters(emale_cogs) +
        geom_seq() +
        geom_bin_label() +
        #geom_feat(size = 5) +
        geom_gene(aes(fill = cluster_id)) +
        scale_fill_brewer("Conserved genes", palette="Set3") +
        #geom_link(color = NA) +
        #scale_fill_distiller(palette = "Spectral") +
        theme()


    emale_cogs %<>% dplyr::mutate(
        cluster_label = paste0(cluster_id, " (", cluster_n, ")"),
        cluster_label = fct_lump_min(cluster_label, 5, other_level = "rare"),
        cluster_label = fct_lump_min(cluster_label, 15, other_level = "medium"),
        cluster_label = fct_relevel(cluster_label, "rare", after=Inf))


    #

    gggenomes(
        genes = emale_genes, seqs = emale_seqs, links = emale_ava,
        feats = list(emale_tirs, ngaros=emale_ngaros, gc=emale_gc)) |>
        add_sublinks(emale_prot_ava) |>
        sync() + # synchronize genome directions based on links
        geom_feat(position="identity", size=6) +
        geom_seq() +
        geom_link(data=links(2)) +
        geom_bin_label() +
        geom_gene(aes(fill=name)) +
        geom_gene_tag(aes(label=name), nudge_y=0.1, check_overlap = TRUE) +
        geom_feat(data=feats(ngaros), alpha=.3, size=10, position="identity") +
        geom_feat_note(aes(label="Ngaro-transposon"), data=feats(ngaros),
                       nudge_y=.1, vjust=0) +
        geom_wiggle(aes(z=score, linetype="GC-content"), feats(gc),
                    fill="lavenderblush4", position=position_nudge(y=-.2), height = .2) +
        scale_fill_brewer("Genes", palette="Dark2", na.value="cornsilk3")
}
