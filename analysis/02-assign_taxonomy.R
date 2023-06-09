#' This scripts assigns taxonomy to aligned 16S of rhizobia strains

library(tidyverse)
library(rRDPData)
library(rRDP)
library(Biostrings)
source(here::here("analysis/00-metadata.R"))

# Read the mapping file
isolates_ID <- read_csv(paste0(folder_data, "raw/rhizobia/02-sequencing/isolates_for_seq.csv"), col_types = cols()) %>%
    mutate(ID = as.numeric(`Sample Name`)) %>%
    select(ExpID, ID) %>%
    mutate(ExpID = str_replace(ExpID, " ", ""))

# Read 16S sequence
isolates_16S <- read_csv(paste0(folder_data, "temp/01-isolates_16S.csv"), col_types = cols()) %>%
    right_join(isolates_ID, by = "ID") %>%
    select(ExpID, ID, Sequence) %>%
    filter(!is.na(Sequence))

# Make DNA string set object
isolates_seq_set <- DNAStringSet(isolates_16S$Sequence)
names(isolates_seq_set) <- isolates_16S$ExpID

# Use rdp for classification (this needs package rRDPData)
pred <- predict(rdp(), isolates_seq_set, confidence = 0)
conf_score <- attr(pred, "confidence") %>% as.data.frame()
colnames(pred) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
colnames(conf_score) <- paste0(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), "Score")
pred$ExpID <- rownames(pred)
conf_score$ExpID <- rownames(pred)

# Join the predicted taxonomy and isolates information
isolates_RDP <- left_join(isolates_16S, pred, by = "ExpID") %>% left_join(conf_score, "ExpID")
write_csv(isolates_RDP, paste0(folder_data, "temp/02-isolates_RDP.csv"))

# Select Rhizobiaceae Family ----
isolates_rhizo <- isolates_RDP %>%
    #filter(Genus == "Ensifer") %>%
    filter(Family == "Rhizobiaceae") %>%
    mutate(Site = str_sub(ExpID, 1, 1)) %>%
    filter(Site %in% c("H", "L"))

isolates_rhizo %>%
    tabyl(Site) # 8 H and 11 L

write_csv(isolates_rhizo, paste0(folder_data, "temp/02-isolates_rhizo.csv"))

# Check the used rhizobia ----
isolates_RDP %>%
    left_join(isolates_ID) %>%
    filter(Owner == "CYC", Genus == "Ensifer")

isolates_RDP_inocula <- isolates_RDP %>%
    left_join(isolates_ID, by = c("ExpID", "ID")) %>%
    filter(Owner == "CYC", Genus == "Ensifer") %>%
    filter(ExpID %in% c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2")) %>%
    select(Owner, Site, ExpID, ID, everything(), Sequence)

write_csv(isolates_RDP_inocula, paste0(folder_data, "temp/02-isolates_RDP_inocula.csv"))
