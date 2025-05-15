#' This script assigns the taxonomy

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
names_blast <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# 1. Aggregate the 16s blast results ----
# Read the reference
ref_16s_seq <- Biostrings::readDNAStringSet("~/bioinformatics/16s/refseq_16s.fasta")
ref_16s <- tibble(accession = str_remove(names(ref_16s_seq), "\\s.+"), scomment = str_remove(names(ref_16s_seq), "^[A-Z]+_[\\d]+.\\d\\s"))

# Read the blast results
list_b_16s <- rep(list(NA), nrow(isolates))
for (i in 1:nrow(isolates)) {
    list_b_16s[[i]] <- read_table(paste0(folder_data, "genomics/taxonomy/", isolates$genome_id[i],"/16s/blast_16s.txt"), col_names = names_blast) %>%
        mutate(genome_id = isolates$genome_id[i])
}

blast_16s <- bind_rows(list_b_16s) %>%
    group_by(qseqid) %>%
    arrange(desc(bitscore))

top_16s <- blast_16s %>%
    # Find the top hit
    slice(1) %>%
    mutate(accession = sseqid) %>%
    left_join(ref_16s) %>%
    select(genome_id, qseqid, scomment, everything()) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id)

write_csv(blast_16s, paste0(folder_data, "genomics_analysis/taxonomy/blast_16s.csv"))
write_csv(top_16s, paste0(folder_data, "genomics_analysis/taxonomy/top_16s.csv"))

# 2. Aggregate the genome blast results ----
ref_genome_seq <- Biostrings::readDNAStringSet(paste0(folder_data, "genomics/blast_db/genomes.fasta"))
ref_genome <- tibble(accession = str_remove(names(ref_genome_seq), "\\s.+"), scomment = str_remove(names(ref_genome_seq), "^[A-Z]+_[A-Z\\d]+.\\d\\s"))
ref_genome <- ref_genome %>%
    rowwise() %>%
    mutate(species = str_split(scomment, " ")[[1]] %>% `[`(2)) %>%
    mutate(strain = str_remove(scomment, "Sinorhizobium \\w+ |Ensifer \\w+ ") %>% str_remove("strain ") %>% str_remove(",")) %>%
    mutate(strain = str_split(strain, " ")[[1]] %>% `[`(1)) %>%
    mutate(replicon = case_when(
        str_detect(scomment, "chromosome") ~ "chromosome",
        str_detect(scomment, "ctg") ~ str_extract(scomment, "ctg\\d+"),
        str_detect(scomment, "plasmid") ~ str_extract(scomment, "plasmid [A-Z|0-9|a-z|_]+")
    ))

# Manually correct the comments
ref_genome$replicon[ref_genome$accession == "NC_003047.1"] <- "chromosome" # NC_003047.1 is 1021 chromsome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006965.1/
ref_genome$replicon[ref_genome$accession == "NC_009636.1"] <- "chromosome" # NC_009636.1 is WSM419 chromosome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000017145.1/
ref_genome$replicon[ref_genome$accession == "NC_012587.1"] <- "chromosome" # NC_012587.1 is NGR234 chromosome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000018545.1/

# Read the blast results
list_b_genome <- rep(list(NA), nrow(isolates))
for (i in 1:nrow(isolates)) {
    list_b_genome[[i]] <- read_table(paste0(folder_data, "genomics/taxonomy/", isolates$genome_id[i],"/blast_genome/blast_genome.txt"), col_names = names_blast) %>%
        mutate(genome_id = isolates$genome_id[i])
}
list_b_genome <- list_b_genome[!is.na(list_b_genome)]

blast_genome <- bind_rows(list_b_genome) %>%
    #filter(pident > 90, bitscore > 10000) %>%
    mutate(accession = sseqid) %>%
    left_join(ref_genome) %>%
    ungroup() %>%
    select(genome_id, species, strain, replicon, qseqid, scomment, sseqid, pident, length, bitscore) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id, desc(bitscore)) %>%
    # Remove those that are symbiotic species but blast to nonsymbiotic plasmids
    #filter(!(!genome_id %in% c("g2", "g3", "g15") & species %in% c("canadensis", "adhaerens"))) %>%
    rename(contig_id = qseqid) %>%
    mutate(contig_id = paste0(genome_id, "_", contig_id)) %>%
    select(genome_id, contig_id, species, strain, replicon, pident, length, bitscore)

top_genome <- blast_genome %>%
    group_by(genome_id, contig_id) %>%
    arrange(desc(bitscore)) %>%
    # Find the top hit
    slice(1)

write_csv(ref_genome, paste0(folder_data, "genomics_analysis/taxonomy/ref_genome.csv"))
write_csv(blast_genome, paste0(folder_data, "genomics_analysis/taxonomy/blast_genome.csv"))
write_csv(top_genome, paste0(folder_data, "genomics_analysis/taxonomy/top_genome.csv"))

