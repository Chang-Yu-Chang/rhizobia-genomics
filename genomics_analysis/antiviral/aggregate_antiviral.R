#' This script aggregate the antiviral data

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

# defense finder
list_def <- list()
for (i in 1:nrow(isolates)) {
    list_def[[isolates$genome_id[i]]] <- read_tsv(paste0(folder_data, "genomics/antiviral/", isolates$genome_id[i], "/defense_finder/", isolates$genome_id[i], "_defense_finder_genes.tsv"), show_col_types = F)
    #read_tsv(paste0(folder_data, "genomics/antiviral/", isolates$genome_id[i], "/defense_finder/", isolates$genome_id[i], "_defense_finder_systems.tsv"), show_col_types = F)
}
defensefinder <- bind_rows(list_def, .id = "genome_id")
write_csv(defensefinder, paste0(folder_data, "genomics_analysis/antiviral/defensefinder.csv"))


# padloc
list_padloc <- list()
for (i in 1:nrow(isolates)) list_padloc[[isolates$genome_id[i]]] <- read_csv(paste0(folder_data, "genomics/antiviral/", isolates$genome_id[i], "/padloc/", isolates$genome_id[i], "_padloc.csv"), show_col_types = F)
padloc <- bind_rows(list_padloc, .id = "genome_id")
write_csv(padloc, paste0(folder_data, "genomics_analysis/antiviral/padloc.csv"))

