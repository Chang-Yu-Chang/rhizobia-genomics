#' This script computes Fst for a set of genes

renv::load()
library(tidyverse)
library(apex) # for reading multiple fasta
library(poppr) # for processing fasta
library(mmod) # for computing Fst
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))

compute_by_gradient <- function (gidi) {
    #' Compute the fst by gradient. gidi is a subset from gid
    #gidi <- gid1
    n_alleles <- nAll(gidi)
    loci_to_keep <- names(n_alleles)[n_alleles >= 2]
    ind <- str_remove(colnames(gidi@tab), "\\.\\w+") %in% loci_to_keep
    gidi <- gidi[,ind] # remove monophoric loci

    dstat <- diff_stats(gidi) # mmod function to differentiation stats
    tb_snp_locus <- bind_cols(tibble(snp_id = rownames(dstat$per.locus)), as_tibble(dstat$per.locus))
    tb_snp_global <- dstat$global
    return(list(tb_snp_locus = tb_snp_locus, tb_snp_global = tb_snp_global))

}


# Set 1. 988 core genes shared by all 36 genomes ----
list_sccg <- gpa$gene[which(apply(gpa[,-1], 1, sum) == 36)]
tb_sccg_fst <- tibble(gene = list_sccg, fst = NA)
folder_set <- paste0(folder_data, "genomics_analysis/fst/set1/")
if (!dir.exists(folder_set)) dir.create(path = folder_set, showWarnings = T)


for (i in 1:length(list_sccg)) {
    cat("\n",i , list_sccg[i])
    # Read alignment fasta
    aln_file <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences/", list_sccg[i], ".aln.fas")
    fasta_data <- read.multiFASTA(aln_file)
    #pegas::nuc.div(fasta_data@dna$group_19086.aln.fas)
    # Correct index names
    gid <- multidna2genind(fasta_data)
    indNames(gid) <- str_remove(indNames(gid), ";\\w+") %>% str_remove("_R_")
    if (length(indNames(gid)) > nrow(genomes)) {cat("\ngene ", i, ": ", list_sccg[i], " is not single copy"); next}
    # Remove the paralog. Only chose the first one -> single copy
    #gid <- gid[!duplicated(indNames(gid)), ]

    # Assign populations
    isolates_pop <- tibble(genome_id = indNames(gid)) %>% left_join(isolates, by = join_by(genome_id))
    strata(gid) <- isolates_pop
    setPop(gid) <- ~site_group
    #slot(gid, "loc.names")
    # diff_stats(gid)
    # diff_test(gid)

    # Subset by pop
    cat("\televation")
    gid1 <- popsub(gid, c("high elevation", "low elevation"))
    tb_snp1 <- compute_by_gradient(gid1)
    cat("\turbanization")
    gid2 <- popsub(gid, c("suburban", "urban"))
    tb_snp2 <- compute_by_gradient(gid2)

    tb_snp <- bind_rows(
        bind_cols(tibble(gradient = "elevation"), tb_snp1$tb_snp_locus),
        bind_cols(tibble(gradient = "urbanization"), tb_snp2$tb_snp_locus)
    ) %>%
        mutate(across(3:7, ~ round(.x, 4)))

    write_csv(tb_snp, paste0(folder_set, list_sccg[i], "-snp.csv"))

    # All 4 populations
    tb_sccg_fst$fst[i] <- list(tibble(metric = names(tb_snp1$tb_snp_global), elev = tb_snp1$tb_snp_global, urba = tb_snp2$tb_snp_global))

    # All pairwise pop across gradients
    tb_gst <- pairwise_Gst_Nei(gid) %>%
        as.matrix() %>% as_tibble() %>%
        mutate(pop1 = colnames(.)) %>%
        pivot_longer(-pop1, names_to = "pop2", values_to = "Gst_est") %>%
        filter(pop1 > pop2) %>%
        mutate(Gst_est = round(Gst_est, 4))
    write_csv(tb_gst, paste0(folder_set, list_sccg[i], "-pop.csv"))
}

# Is single copy?
tb_sccg_fst$singlecopy <- NA

for (i in 1:length(list_sccg)) {
    cat("\n",i , list_sccg[i])
    # Read alignment fasta
    aln_file <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences/", list_sccg[i], ".aln.fas")
    fasta_data <- read.multiFASTA(aln_file)
    # Correct index names
    gid <- multidna2genind(fasta_data)
    indNames(gid) <- str_remove(indNames(gid), ";\\w+") %>% str_remove("_R_")
    if (length(indNames(gid)) > nrow(genomes)) {
        tb_sccg_fst$singlecopy[i] <- F
        cat("\ngene ", i, ": ", list_sccg[i], " is not single copy")
    } else if (length(indNames(gid)) == nrow(genomes)){
        tb_sccg_fst$singlecopy[i] <- T
    }
}

sccg_fst <- tb_sccg_fst %>% unnest(cols = fst)

write_csv(sccg_fst, paste0(folder_data, "genomics_analysis/fst/set1_fst.csv"))



# Set 2. core genes shared by 32 symbiotic strains ----
gpa2 <- gpa %>% select(-g2, -g3, -g15, -g42) # remove g42 for its bad quality
list_sccg <- gpa2$gene[which(apply(gpa2[,-1], 1, sum) == 32)]
tb_sccg_fst <- tibble(gene = list_sccg, fst = NA)
folder_set <- paste0(folder_data, "genomics_analysis/fst/set2/")
if (!dir.exists(folder_set)) dir.create(path = folder_set, showWarnings = T)


for (i in 1:length(list_sccg)) {
    cat("\n",i , list_sccg[i])
    # Read alignment fasta
    aln_file <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences/", list_sccg[i], ".aln.fas")
    fasta_data <- read.multiFASTA(aln_file)
    # Correct index names
    gid <- multidna2genind(fasta_data)
    indNames(gid) <- str_remove(indNames(gid), ";\\w+") %>% str_remove("_R_")
    if (length(indNames(gid)) > nrow(genomes)) {cat("\ngene ", i, ": ", list_sccg[i], " is not single copy"); next}
    # Remove the paralog. Only chose the first one -> single copy
    #gid <- gid[!duplicated(indNames(gid)), ]

    # Assign populations
    isolates_pop <- tibble(genome_id = indNames(gid)) %>% left_join(isolates, by = join_by(genome_id))
    strata(gid) <- isolates_pop
    setPop(gid) <- ~site_group
    #slot(gid, "loc.names")
    # diff_stats(gid)
    # diff_test(gid)

    # Subset by pop
    cat("\televation")
    gid1 <- popsub(gid, c("high elevation", "low elevation"))
    tb_snp1 <- compute_by_gradient(gid1)
    cat("\turbanization")
    gid2 <- popsub(gid, c("suburban", "urban"))
    tb_snp2 <- compute_by_gradient(gid2)

    tb_snp <- bind_rows(
        bind_cols(tibble(gradient = "elevation"), tb_snp1$tb_snp_locus),
        bind_cols(tibble(gradient = "urbanization"), tb_snp2$tb_snp_locus)
    ) %>%
        mutate(across(3:7, ~ round(.x, 4)))

    write_csv(tb_snp, paste0(folder_set, list_sccg[i], "-snp.csv"))

    # All 4 populations
    tb_sccg_fst$fst[i] <- list(tibble(metric = names(tb_snp1$tb_snp_global), elev = tb_snp1$tb_snp_global, urba = tb_snp2$tb_snp_global))

    # All pairwise pop across gradients
    tb_gst <- pairwise_Gst_Nei(gid) %>%
        as.matrix() %>% as_tibble() %>%
        mutate(pop1 = colnames(.)) %>%
        pivot_longer(-pop1, names_to = "pop2", values_to = "Gst_est") %>%
        filter(pop1 > pop2) %>%
        mutate(Gst_est = round(Gst_est, 4))
    write_csv(tb_gst, paste0(folder_set, list_sccg[i], "-pop.csv"))
}

# Is single copy ?
tb_sccg_fst$singlecopy <- NA

for (i in 1:length(list_sccg)) {
    cat("\n",i , list_sccg[i])
    # Read alignment fasta
    aln_file <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences/", list_sccg[i], ".aln.fas")
    fasta_data <- read.multiFASTA(aln_file)
    # Correct index names
    gid <- multidna2genind(fasta_data)
    indNames(gid) <- str_remove(indNames(gid), ";\\w+") %>% str_remove("_R_")
    if (length(indNames(gid)) > nrow(genomes)) {
        tb_sccg_fst$singlecopy[i] <- F
        cat("\ngene ", i, ": ", list_sccg[i], " is not single copy")
    } else if (length(indNames(gid)) == nrow(genomes)){
        tb_sccg_fst$singlecopy[i] <- T
    }
}

sccg_fst <- tb_sccg_fst %>% unnest(cols = fst)

write_csv(sccg_fst, paste0(folder_data, "genomics_analysis/fst/set2_fst.csv"))







if (F) {

#
aln_file <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences/", list_sccg[i], ".aln.fas")
fasta_data <- read.multiFASTA(aln_file)
gid1 <- popsub(gid, c("high elevation", "low elevation"))
dstat_snps <- diff_stats(gid1)$per.locus
dstat_gene <- diff_stats(gid1)$global

dtest_snps <- diff_test(gid1, sim = T, nreps = 1000)

mean(dtest_snps)
bs <- chao_bootstrap(gid1, nreps = 10)

aln_dist <- dist.multidna(fasta_data, pool = TRUE)
aln_dist
amova(aln_dist ~ population, data = strata(gid), nperm = 100)

# Plot the gene
dstat_snps %>%
    as_tibble() %>%
    mutate(snp_id = as.numeric(rownames(dstat_snps))) %>%
    mutate(p_value = dtest_snps) %>%
    mutate(signif = p_value < 0.05) %>%
    ggplot() +
    geom_point(aes(x = snp_id, y = Gst, color = signif)) +
    scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "red")) +
    theme_bw() +
    theme() +
    guides() +
    labs()


# xx <- pairwise_Gst_Nei(gid)
# mean(xx)


















}
