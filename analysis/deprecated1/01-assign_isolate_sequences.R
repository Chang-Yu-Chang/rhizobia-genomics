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
    REGEX_SuffixReverse = "R.ab1$", M1TrimmingCutoff = 1e-3
)

#launchApp(alignment1)
#read.abif(paste0(folder_data, "raw/sanger/1st Round/Wood06_68-16S-rRNA-seqF.ab1"))

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
    REGEX_SuffixReverse = "R.ab1$", M1TrimmingCutoff = 1e-3
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

