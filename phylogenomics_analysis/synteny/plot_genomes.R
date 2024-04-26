#' This script plots the genomes

renv::load()
library(tidyverse)
library(gggenomes)
source(here::here("metadata.R"))

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



if (FALSE) {
data(package="gggenomes")

emale_seqs
# emale_seqs <- emale_seqs %>%
#   extract(seq_desc, into = c("emale_type", "is_typespecies"), "=(\\S+) \\S+=(\\S+)", remove=F, convert=T) %>%
#   dplyr::arrange(emale_type, length)
emale_seqs[1:6,]

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
#
