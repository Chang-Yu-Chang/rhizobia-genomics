#' This scripts works on the genome data from plasmidsaurus

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# Read the mapping files
isolates_bgs <- read_csv(here::here("protocols/20230830_plasmidsaurus.csv"), col_names = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv")) %>%
    filter(ID %in% isolates_bgs$`X1`)


# Read the plasmidsaurus analysis
list_raw <- c(paste0("Chang_Q5C_results/Chang_Q5C_", c(1:10, 12:17, 19)), paste0("Chang_Q5C_results_repeated/Chang_Q5C_", c(11,18)))

## Does the sourmash predicted taxonomy match the Sanger RDP?
list_sourmash <- rep(list(NA), length(list_raw))
for (i in 1:length(list_raw)) {
    list_sourmash[[i]] <- read_csv(paste0(folder_data, "raw/", list_raw[i], "/analysis/sourmash_species.txt"), show_col_types = F) %>%
        mutate(strain_ID = c(1:10, 12:17, 19, 11, 18)[i])
}

list_sourmash %>%
    bind_rows() %>%
    filter(!is.na(genus)) %>%
    select(strain_ID, everything()) %>%
    group_by(strain_ID) %>%
    slice(1)

isolates_RDP
# Isolate ID 1, 18, 48, 54, 60 are rhizobium, which match the Sanger RDP results


##
