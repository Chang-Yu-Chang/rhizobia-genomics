#' This script analyze the called VCF file
library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
source(here::here("analysis/00-metadata.R"))

# 0. read data ----
#g2_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_2/snps.vcf"))
#g3_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_3/snps.vcf"))
#g4_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_4/snps.vcf"))
# fb_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/freebayes/test.vcf"))
# fb_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/freebayes/var.vcf"))
#snf_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/sniffles/sv.vcf")) # structure variants
core_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/core/core.vcf"))
medicae_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy_medicae/core/core.vcf"))

vcfR::getCHROM(medicae_vcf) %>% table() %>% sum()

# Loaded 3 sequences totalling 6716226 bp.
# Will mask 0 regions totalling 0 bp ~ 0.00%
# 0	Chang_Q5C_2	        snp=34029	del=92	    ins=98	het=325	unaligned=5719454
# 1	Chang_Q5C_3 	    snp=34328	del=116	    ins=117	het=1202	unaligned=5695322
# 2	Chang_Q5C_4     	snp=243405	del=1561	ins=1616	het=3324	unaligned=2656123
# 3	Chang_Q5C_5	        snp=245029	del=1597	ins=1662	het=3332	unaligned=2676579
# 4	Chang_Q5C_6  	    snp=241758	del=1539	ins=1602	het=5163	unaligned=2629813
# 5	Chang_Q5C_8	        snp=242787	del=1655	ins=1666	het=3350	unaligned=2672835
# 6	Chang_Q5C_9	        snp=243677	del=1598	ins=1674	het=3095	unaligned=2635316
# 7	Chang_Q5C_10	    snp=29542	del=2050	ins=2262	het=3755	unaligned=705034
# 8	Chang_Q5C_11	    snp=240940	del=1606	ins=1638	het=5162	unaligned=2638271
# 9	Chang_Q5C_13    	snp=242648	del=1585	ins=1703	het=5021	unaligned=2663305
# 10	Chang_Q5C_15	snp=40148	del=119	    ins=157	    het=1376	unaligned=5595621
# 11	Chang_Q5C_16	snp=244097	del=1544	ins=1658	het=4884	unaligned=2650581
# 12	Chang_Q5C_17	snp=242988	del=1554	ins=1674	het=4756	unaligned=2634310
# 13	Chang_Q5C_19	snp=244143	del=1681	ins=1709	het=3760	unaligned=2554793

# Check the taxonomy
isolates <- read_csv(paste0(folder_data, "temp/42-isolates.csv"), show_col_types = F) %>%
    mutate(genome_name = str_replace(genome_id, "g", "Chang_Q5C_") %>% factor(paste0("Chang_Q5C_", 1:20)))

mash_g_top <- read_csv(paste0(folder_data, "temp/38-mash_g_top.csv"))
mash_g_top %>%
    select(genome_id, identity, query_comment) %>%
    #filter(str_detect(query_comment, "Ensifer|Sinorhizobium")) %>%
    select(-identity) %>%
    group_by(genome_id) %>%
    mutate(mash_hits = paste0("mash", 1:n())) %>%
    ungroup() %>%
    pivot_wider(id_cols = genome_id, names_from = mash_hits, values_from = query_comment) %>%
    left_join(select(isolates, genome_id, strain_site_group)) %>%
    select(site = strain_site_group, genome_id, everything())
# g2, 3, 10, 15 are meliloti
#  the rest should be medicae

# Convert the vcf to a genclone object
#core_genind <- vcfR2genind(core_vcf)
core_genind <- vcfR2genind(medicae_vcf)
# core_genind <- as.genclone(core_genind)
isolates_ensifer <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F) %>%
    filter(genus == "Ensifer") %>%
    mutate(genome_name = str_replace(genome_id, "g", "Chang_Q5C_")) %>%
    select(genome_name, everything()) %>%
    mutate(genome_name = factor(genome_name, rownames(core_genind$tab))) %>%
    arrange(genome_name) %>%
    drop_na()
other(core_genind)$taxa <- isolates_ensifer

# Assign populations
strata(core_genind) <- data.frame(other(core_genind)$taxa)
setPop(core_genind) <- ~ site
popNames(core_genind)
length(names(core_genind$all.names))
names(core_genind$all.names) %>%
    str_split("_") %>%
    sapply(function(x)`[`(x,2))  %>%
    table()


# 1. remove low-quality SNPs in the vcf
#vcf_data[core_vcf$QUAL >= 30, ]

# 1. check snippy results ----
# Remove uninformative loci
core_genind_cut <- informloci(core_genind, cutoff = 2/nInd(core_genind), MAF = 0.01, quiet = FALSE)
popNames(core_genind_cut)
#core_genind
core_genind_cut

# mlg
mlg(core_genind_cut)
mlg.crosspop(core_genind_cut, quiet = F)

# nucleotide diversity
mlg_tab <- mlg.table(core_genind_cut, plot = F)
diversity_stats(mlg_tab)
# Index
# Pop        H G    lambda E.5
# L 2.079442 8 0.8750000   1
# H 1.791759 6 0.8333333   1
# H	 logical whether or not to calculate Shannon's index
# G	 logical whether or not to calculate Stoddart and Taylor's index (aka inverse Simpson's index).
# lambda logical whether or not to calculate Simpson's index
# E5 logical whether or not to calculate Evenness

# Summary
poppr(core_genind)
# Pop  N MLG eMLG SE    H  G lambda E.5   Hexp    Ia rbarD        File
# 1     L  8   8    8  0 2.08  8  0.875   1 0.0247 15479 0.589 core_genind
# 2     H  6   6    6  0 1.79  6  0.833   1 0.0691 17789 0.671 core_genind
# 3 Total 14  14   10  0 2.64 14  0.929   1 0.0857 15509 0.465 core_genind

# 1. PCA using the 14 Ensifer strains ----
# PCA
core_genlight_cut <- dartR::gi2gl(core_genind_cut)
core_pca <- glPca(core_genlight_cut, nf = 2, center = T, scale = T)
eig <- sort(core_pca$eig/sum(core_pca$eig), decreasing = T)

# Merge the tb
isolates_pca <- core_pca$scores %>%
    as_tibble() %>%
    mutate(PC1 = PC1/100, PC2 = PC2/100) %>%
    mutate(genome_name = core_genlight_cut$ind.names) %>%
    left_join(isolates_ensifer)

set.seed(1)
p1 <- isolates_pca %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site), shape = 21, size = 3, stroke = 1.5, position = position_jitter(width = 0.02, height = 0.02)) +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        #legend.background = element_rect(fill = "white", color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = 2),
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(eig[1]*100,1), "%)"),
         y = paste0("PC2 (", round(eig[2]*100,1), "%)"))

# include the mash result
p2 <- isolates_pca %>%
    # THIS IS HARD CODED BY CHECKING THE MASH OUTCOME. COME UP WITH A BETTER WAY
    left_join(tibble(genome_name = paste0("Chang_Q5C_", c(2,3,10,15,4:6, 8,9,11,13,16,17,19)),
                     species_mash = c(rep("meliloti", 4), rep("medicae", 10)))) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = species_mash), shape = 21, size = 3, stroke = 1.5, position = position_jitter(width = 0.02, height = 0.02)) +
    scale_color_npg() +
    theme_classic() +
    theme(
        #legend.background = element_rect(fill = "white", color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = 2),
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(eig[1]*100,1), "%)"),
         y = paste0("PC2 (", round(eig[2]*100,1), "%)"))

p <- plot_grid(p1, p2, nrow = 1, align = "h", scale = 0.95, axis = "tb", labels = c("A", "B")) + paint_white_background()

ggsave(paste0(folder_data, "temp/47-01-snippy_core.png"), width = 8, height = 4.5)

# 2. PCA using the 10 Ensifer strains ----
# PCA
core_genlight_cut <- dartR::gi2gl(core_genind_cut)
core_pca <- glPca(core_genlight_cut, nf = 2, center = T, scale = T)
eig <- sort(core_pca$eig/sum(core_pca$eig), decreasing = T)

# Merge the tb
isolates_pca <- core_pca$scores %>%
    as_tibble() %>%
    mutate(PC1 = PC1/100, PC2 = PC2/100) %>%
    mutate(genome_name = core_genlight_cut$ind.names) %>%
    left_join(isolates_ensifer)

set.seed(1)
p1 <- isolates_pca %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site), shape = 21, size = 3, stroke = 1.5, position = position_jitter(width = 0.02, height = 0.02)) +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        #legend.background = element_rect(fill = "white", color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = 2),
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(eig[1]*100,1), "%)"),
         y = paste0("PC2 (", round(eig[2]*100,1), "%)"))

# include the mash result
p2 <- isolates_pca %>%
    # THIS IS HARD CODED BY CHECKING THE MASH OUTCOME. COME UP WITH A BETTER WAY
    left_join(tibble(genome_name = paste0("Chang_Q5C_", c(2,3,10,15,4:6, 8,9,11,13,16,17,19)),
                     species_mash = c(rep("meliloti", 4), rep("medicae", 10)))) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = species_mash), shape = 21, size = 3, stroke = 1.5, position = position_jitter(width = 0.02, height = 0.02)) +
    scale_color_npg() +
    theme_classic() +
    theme(
        #legend.background = element_rect(fill = "white", color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = 2),
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(eig[1]*100,1), "%)"),
         y = paste0("PC2 (", round(eig[2]*100,1), "%)"))

p <- plot_grid(p1, p2, nrow = 1, align = "h", scale = 0.95, axis = "tb", labels = c("A", "B")) + paint_white_background()

ggsave(paste0(folder_data, "temp/47-02-snippy_medicae.png"), width = 8, height = 4.5)







if (F) {



# rarefaction
# library(vegan)
# mlg_tab <- mlg.table(core_genind, plot = F)
# rarecurve(mlg_tab, ylab="Number of expected MLGs", sample=min(rowSums(mlg_tab)),
#           border = NA, fill = NA, font = 2, cex = 1, col = "blue")




# Check data
# missingno(core_genind) # Detect missing data
# info_table(core_genind, plot = T)
nInd(core_genind) # Number of individuals = 14
locNames(core_genind) # Names of the loci

#
popsub(core_genind, sublist = c())

data("H3N2", package = "adegenet")
strata(H3N2) <- data.frame(other(H3N2)$x)
H3N2
#nInd(H3N2)
other(H3N2)$x %>% as_tibble()
setPop(H3N2) <- ~country
v_na <- popsub(H3N2, sublist = c("USA", "Canada"))
popNames(v_na)
c(NorthAmerica = nInd(v_na), Total = nInd(H3N2))
v_na_minus <- popsub(H3N2, exclude = c("USA", "Canada"))
popNames(v_na_minus)
(nInd(v_na_minus) + nInd(v_na)) == nInd(H3N2)

# nLoc(core_genind)
# core_genind@type
# core@samples <- merge(vcf_data@samples, population_mapping, by = "SampleID", all.x = TRUE)


# 0.1 summary
core_vcf
# ***** Object of Class vcfR *****
# 14 samples
# 3 CHROMs
# 34,156 variants
# Object size: 7.8 Mb
# 0 percent missing data
# *****        *****         *****
core_genind
#' /// GENIND OBJECT /////////
#'
#' // 14 individuals; 34,156 loci; 73,143 alleles; size: 25.4 Mb
#'
#' // Basic content
#' @tab:  14 x 73143 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 1-4)
#' @loc.fac: locus factor for the 73143 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: adegenet::df2genind(X = t(x), sep = sep)

# 0.2 Add ----
vcf_data@samples <- merge(vcf_data@samples, population_mapping, by = "SampleID", all.x = TRUE)

# Allele frequencies
poppr_result <- poppr(core_genind, method = "freq")
poppr_result
# Pop  N MLG eMLG SE    H  G lambda E.5   Hexp    Ia rbarD        File
# 1 Total 14  14   14  0 2.64 14  0.929   1 0.0859 15532 0.465 core_genind

# 1. plot snps ----
snps_wide <- as_tibble(core_genind$tab)
snps <- snps_wide %>%
    mutate(genome_name = paste0("Chang_Q5C_", c(2, 3, 4, 5, 6 , 8, 9, 10, 11, 13, 15, 16, 17, 19))) %>%
    pivot_longer(cols = -genome_name, names_to = "locus_allele", values_to = "value") %>%
    filter(value == 1) %>%
    separate(col = locus_allele, into = c("locus", "allele"), sep = "\\.") %>%
    select(-value) %>%
    mutate(contig = str_sub(locus, 1, 14)) %>%
    mutate(allele = factor(allele))

p <- snps %>%
    ggplot() +
    geom_tile(aes(x = locus, y = genome_name, fill = allele)) +
    #scale_fill_manual(values = c(`1` = "maroon", `0` = "snow")) +
    scale_fill_aaas() +
    facet_grid(.~contig, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_blank()
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/47-01-snps.png"), p, width = 30, height = 10)

# Numbers
snps %>%
    distinct(locus, .keep_all = T) %>%
    group_by(contig) %>%
    count()
# contig             n
# <chr>          <int>
#     1 c_000000000001 30963
# 2 c_000000000002    69
# 3 c_000000000003  3124


# pi_values <- pi_result$pi
# The pi_values variable now contains the









# ref_dna <- ape::read.dna(paste0(folder_data, "temp/plasmidsaurus/usda1106/genome.fasta"), format = "fasta")
# ref_gff <-  read.table(paste0(folder_data, "temp/plasmidsaurus/usda1106/genome.gff"), sep="\t", quote="")
#core_chrom <- create.chromR(vcf = core_vcf, seq = ref_dna, ann = ref_gff)
}
