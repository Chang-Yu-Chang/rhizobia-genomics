#' This script plots the tree based on the 16S


library(tidyverse)
#library(sangeranalyseR)
library(msa)
library(tools)
source(here::here("analysis/00-metadata.R"))

isolates_RDP <- read_csv(paste0(folder_data, "temp/01-isolates_RDP.csv"), show_col_types = F)
gc_plate <- read_csv(paste0(folder_data, "raw/growth_curve/gc_plate.csv"))
list_strains <- unique(gc_plate$strain)[-20]

list_sequences <- isolates_RDP %>%
    filter(ExpID %in% list_strains) %>%
    pull(Sequence) %>%
    DNAStringSet()

aln <- msa(list_sequences, "ClustalW")
print(aln, show="complete")
# msaPrettyPrint(aln, output="tex", showNames="none", showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
# texi2dvi("aln.tex", clean=TRUE)

