#' This script plot the Fst per gene and per snp

library(tidyverse)
library(cowplot)
library(ggsci)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))

#
set_name
gene_wide_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/elev_med/gene_wide_fst.csv"))
hist(gene_wide_fst$Gst_est)

range(gene_wide_fst$Gst_est, na.rm = T)

gene_wide_fst %>%
    filter(Gst_est >= 0.8)

per_locus_fst %>%
    filter(fst.Gst == 1) %>%
    view


# Plot tree
which(list_sccg$gene == "yhdY_3~~~yhdY_2")
i=4643
list_sccg$gene[i]
tr <- read.tree(paste0(folder_data, "phylogenomics_analysis/trees/elev_med/seq_core/", list_sccg$gene[i], "/", list_sccg$gene[i], ".contree"))
plot(tr)
