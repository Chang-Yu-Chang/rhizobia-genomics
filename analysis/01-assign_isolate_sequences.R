#' This scripts reads and aligns raw sanger sequences from Genewiz

if (FALSE) {
    install.packages("BiocManager") # Package for managing and installing Bioconductor packages
    BiocManager::install("sangeranalyseR")
}

library(tidyverse)
library(sangeranalyseR)
folder_data <- here::here("data/")
source(here::here("analysis/00-metadata.R"))


# 1. Align the F and R reads ----
alignment1 <- sangeranalyseR::SangerAlignment(
    ABIF_Directory = paste0(folder_data, "raw/sanger/1st Round/"),
    REGEX_SuffixForward = "F.ab1$",
    REGEX_SuffixReverse = "R.ab1$"
)

#' This generates three merged fasta files:
#' 1. Sanger_all_trimmed_reads.fa
#' 2. Sanger_contigs_alignment.fa
#' 3. Sanger_contigs_unalignment.fa

sangeranalyseR::writeFasta(
    alignment1,
    outputDir = paste0(folder_data, "raw/sanger/1st Round"),
    compress = FALSE,
    compression_level = NA,
    selection = "all"
)

# Repeat it for the 2rd round
alignment2 <- sangeranalyseR::SangerAlignment(
    ABIF_Directory = paste0(folder_data, "raw/sanger/2nd Round/"),
    REGEX_SuffixForward = "F.ab1$",
    REGEX_SuffixReverse = "R.ab1$"
)
sangeranalyseR::writeFasta(
    alignment2,
    outputDir = paste0(folder_data, "raw/sanger/2nd Round"),
    compress = FALSE,
    compression_level = NA,
    selection = "all"
)



#
consensus1 <- seqinr::read.fasta(paste0(folder_data, "raw/sanger/1st Round/Sanger_contigs_alignment.fa"), seqtype = "DNA")
consensus2 <- seqinr::read.fasta(paste0(folder_data, "raw/sanger/2nd Round/Sanger_contigs_alignment.fa"), seqtype = "DNA")
clean_concensus <- function(x) {
    contig <- x %>% as.character()
    contig[contig=="a"] <- "A"
    contig[contig=="t"] <- "T"
    contig[contig=="g"] <- "G"
    contig[contig=="c"] <- "C"
    contig <- paste(contig, collapse = "")
    contig <- gsub("-", "", contig)
    return(data.frame(Sequence = contig))
}
consensus1 <- consensus1 %>% lapply(clean_concensus)
consensus2 <- consensus2 %>% lapply(clean_concensus)

df_seq1 <- bind_rows(consensus1, .id = "temp") %>%
    tidyr::separate(col = "temp", into = c("ID", "t1", "t2", "t3"), sep = "-") %>%
    dplyr::select(ID, Sequence) %>%
    dplyr::mutate(ID = str_replace(ID, "Wood06_", "") %>% as.numeric() %>% factor(levels = 1:100)) %>%
    dplyr::mutate(ConsensusLength = nchar(Sequence)) %>%
    as_tibble() %>%
    arrange(ID)

df_seq2 <- bind_rows(consensus2, .id = "temp") %>%
    tidyr::separate(col = "temp", into = c("ID", "t1", "t2", "t3"), sep = "-") %>%
    dplyr::select(ID, Sequence) %>%
    dplyr::mutate(ID = str_replace(ID, "Wood06_", "") %>% as.numeric() %>% factor(levels = 1:100)) %>%
    dplyr::mutate(ConsensusLength = nchar(Sequence)) %>%
    as_tibble() %>%
    arrange(ID)

df_seq <- filter(df_seq1, !(df_seq1$ID %in% df_seq2$ID)) %>%
    bind_rows(df_seq2) %>%
    arrange(ID) %>%
    distinct(ID, .keep_all = T)

write_csv(df_seq, paste0(folder_data, "temp/01-isolates_16S.csv"))

# 2. RDP ----
#' This scripts reads isolate 16S sequences and assigns taxonomy using RDP
if (FALSE) {
    BiocManager::install("rRDPData")
    BiocManager::install("rRDP")
}

library(tidyverse)
library(rRDPData)
library(rRDP)
library(Biostrings)
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
write_csv(isolates_RDP, paste0(folder_data, "temp/01-isolates_RDP.csv"))


# 3. Plot ----
library(cowplot)
isolates_RDP <- read_csv(paste0(folder_data, "temp/01-isolates_RDP.csv"), show_col_types = F) %>%
    #arrange(Genus) %>%
    # mutate(Genus = factor(Genus, levels = rev(unique(isolates_RDP$Genus)))) %>%
    # mutate(Family = factor(Family, levels = rev(unique(isolates_RDP$Family)))) %>%
    filter(!is.na(Family))

isolates_RDP %>%
    filter(ExpID %in% c("H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1", "L4M2R2", "L4M3R3", "L4M4R1")) %>%
    view
#
isolates_rhizo <- isolates_RDP %>%
    #filter(Genus == "Ensifer") %>%
    filter(Family == "Rhizobiaceae") %>%
    mutate(Site = str_sub(ExpID, 1, 1)) %>%
    filter(Site %in% c("H", "L"))

isolates_rhizo %>%
    tabyl(Site) # 8 H and 11 L

write_csv(isolates_rhizo, paste0(folder_data, "temp/01-isolates_rhizo.csv"))


#
p1 <- isolates_RDP %>%
    ggplot() +
    geom_bar(aes(y = Family, fill = Family), color = 1) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(isolates_RDP)), vjust = 2, hjust = 2) +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = grey(0.8), linetype = 1))
    #guides(fill = "none")

p2 <- isolates_RDP %>%
    ggplot() +
    geom_bar(aes(y = Genus, fill = Family), color = 1) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(isolates_RDP)), vjust = 2, hjust = 2) +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = grey(0.8), linetype = 1))


p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")

ggsave(here::here("plots/01-isolate_taxonomy.png"), plot = p, width = 8, height = 7)

# Isolate ID distribution
isolates_ID <- isolates_ID %>%
    mutate(Site = ExpID %>% str_sub(1, 3) %>% str_replace("-", "") %>% str_replace("M", "") %>%
               str_replace("fp1|fp2", "fp") %>% str_replace("gp1", "gp")) %>%
    mutate(Owner = ifelse(str_sub(Site, 1, 1) %in% c("L", "H"), "CYC", "TPP"))

site_id <- unique(isolates_ID$Site)
family_id <- c("Rhizobiaceae", "Pseudomonadaceae", "Enterobacteriaceae", "Others")


p <- isolates_RDP %>%
    mutate(Family = as.character(Family)) %>%
    mutate(Family = ifelse(Family %in% c("Rhizobiaceae", "Pseudomonadaceae", "Enterobacteriaceae"), Family, "Others")) %>%
    left_join(isolates_ID) %>%
    mutate(Site = factor(Site, site_id)) %>%
    mutate(Family = factor(Family, family_id)) %>%
    group_by(Owner, Site, Family) %>%
    summarize(Count = n()) %>%
    ggplot() +
    geom_col(aes(x = Site, y = Count, fill = Family), color = 1) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(breaks = 1:10) +
    theme_classic()
ggsave(here::here("plots/01-family_count.png"), plot = p, width = 6, height = 3)

# Check the used rhizobia ----
isolates_RDP %>%
    left_join(isolates_ID) %>%
    filter(Owner == "CYC", Genus == "Ensifer")

isolates_RDP_inocula <- isolates_RDP %>%
    left_join(isolates_ID, by = c("ExpID", "ID")) %>%
    filter(Owner == "CYC", Genus == "Ensifer") %>%
    filter(ExpID %in% c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2")) %>%
    select(Owner, Site, ExpID, ID, everything(), Sequence)

write_csv(isolates_RDP_inocula, paste0(folder_data, "temp/01-isolates_RDP_inocula.csv"))


# expid_used <- c("H1M2R3", "H3M1R1", "H4M1R1", "L1M3R2", "L2M2R1", "L3M1R1") # These are not all rhizboia
# isolates_RDP %>%
#     filter(ExpID %in% expid_used)














