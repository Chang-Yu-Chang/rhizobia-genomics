#' This script plot the annotation

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

ann_bak <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_1/05-bakta/consensus.tsv"), skip = 2, show_col_types = F)
ann_pro <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_1/10-prokka/annotated.tsv"), show_col_types = F)


ann_pro %>%
    filter(str_detect(product, "plasmid"))

"prokka does not tell you which contig it's from..."

# library(GenomicRanges)
# library(ggbio)
# data("CRC", package = "biovizBase")
# head(hg19sub)
# ggbio() +
#     circle(hg19sub, geom = "ideo", fill = "gray70") + #Ideogram
#     circle(hg19sub, geom = "scale", size = 2) + #Scale
#     circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3) # label
#
# sequences <- Biostrings::readDNAStringSet(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_1/04-medaka/consensus.fasta"))
# gr <- GRanges(seqnames = names(sequences), ranges = IRanges(start = 1, end = width(sequences)))
#
# ggbio() +
#     circle(gr, geom = "ideo", fill = "gray70") +
#     circle(gr, geom = "scale", size = 2) + #Scale
#     circle(gr, geom = "text", aes(label = seqnames), vjust = 0, size = 3) # label
#gff <- read.gff(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_2/10-prokka/annotated.gff"))
#gff <- read.gff(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_2/05-bakta/consensus.gff3"))
#
