#' This script analyze the called VCF file
library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
source(here::here("analysis/00-metadata.R"))

# 0. read data ----
g2_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_2/snps.vcf"))
g3_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_3/snps.vcf"))
g4_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_4/snps.vcf"))
g5_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_5/snps.vcf"))
g6_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/Chang_Q5C_6/snps.vcf"))
core_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/core/core.vcf"))
#core_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/08-vcf/core.vcf"))
#head(is.polymorphic(core_vcf, na.omit = TRUE))

# 1. remove low-quality SNPs in the vcf
vcfR::me(g2_vcf, g3_vcf)
vcf_data[core_vcf$QUAL >= 30, ]

# Convert the vcf to a genclone object
core_genind <- vcfR2genind(core_vcf)
# core_genind <- as.genclone(core_genind)
isolates_ensifer <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F) %>%
    filter(genus == "Ensifer") %>%
    mutate(genome_name = str_replace(genome_id, "g", "Chang_Q5C_")) %>%
    select(genome_name, everything()) %>%
    mutate(genome_name = factor(genome_name, rownames(core_genind$tab))) %>%
    arrange(genome_name)
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






# 1. Multi-locus genotypes (MLGs) ----
# Removing uninformative loci
core_genind_cut <- informloci(core_genind, cutoff = 2/nInd(core_genind), MAF = 0.01, quiet = FALSE)
popNames(core_genind_cut)

#
mlg(core_genind_cut)
mlg.crosspop(core_genind_cut, quiet = F)

# pi
mlg_tab <- mlg.table(core_genind_cut, plot = F)
diversity_stats(mlg_tab)


#
library(vegan)
mlg_tab <- mlg.table(core_genind_cut, plot = F)
rarecurve(mlg_tab, ylab="Number of expected MLGs", sample=min(rowSums(mlg_tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")


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
snps_wide %>%
    select(starts_with(c("c_000000000001_9497", "c_000000000001_9521"))) %>%
    mutate(genome_name = paste0("Chang_Q5C_", c(2, 3, 4, 5, 6 , 8, 9, 10, 11, 13, 15, 16, 17, 19))) %>%
    pivot_longer(cols = -genome_name, names_to = "locus_allele", values_to = "value") %>%
    filter(value == 1) %>%
    separate(col = locus_allele, into = c("locus", "allele"), sep = "\\.") %>%
    select(-value)


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









