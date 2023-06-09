#' This script does multi sequence alignment and build trees

library(tidyverse)
library(msa)
library(phangorn)
source(here::here("analysis/00-metadata.R"))

isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F)
gc_plate <- read_csv(paste0(folder_data, "raw/growth_curve/gc_plate.csv"), show_col_types = F)
list_strains <- unique(gc_plate$strain)[-20]
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")

# Subset the 19 strains
isolates_meta <- isolates_RDP %>%
    filter(ExpID %in% list_strains) %>%
    mutate(Site = str_sub(ExpID, 1,1)) %>%
    mutate(label = as.character(1:n()))
write_csv(isolates_meta, paste0(folder_data, "temp/03-isolates_meta.csv"))

# Multiple sequence alignment
list_sequences <- isolates_meta %>%
    pull(Sequence) %>%
    DNAStringSet()
sapply(list_sequences, length) %>% range() # consensus length [410, 1555]
aln <- msa(list_sequences, "ClustalW")


# Save tree
rhizo <- as.phyDat(aln)
dm <- dist.ml(rhizo) # create distance map
tree_NJ <- NJ(dm) # create tree
save(tree_NJ, file = paste0(folder_data, "temp/03-tree_NJ.Rdata"))
