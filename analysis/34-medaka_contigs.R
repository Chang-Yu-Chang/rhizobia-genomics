#' This script is to have a overview on the raw read statistics

library(tidyverse)
library(cowplot)
library(seqinr)
source(here::here("analysis/00-metadata.R"))

# 0. read the raw read data
list_g <- paste0(rep("Chang_Q5C_results", 19), "/Chang_Q5C_", 1:19, "/")
list_g[11] <- "Chang_Q5C_results_repeated/Chang_Q5C_11/"
list_g[18] <- "Chang_Q5C_results_repeated/Chang_Q5C_18/"
list_reads <- rep(list(NA), 19)

i=1
fa1 <- read.fasta(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "04-medaka/consensus.fasta"))

fa1$contig_1

"DRAW The CONTIG and ACCUMULATIVE SIZE"
