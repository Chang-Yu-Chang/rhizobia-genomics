#' This script assigns the taxonomy

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
names_blast <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

if (F) {

# 1. Aggregate sourmash ----
list_sm <- rep(list(NA), nrow(isolates))
for (i in 1:length(list_sm)) {
    tb <- read_csv(paste0(folder_genomics, "taxonomy/", isolates$genome_id[i], "/sourmash/gathered.csv"))
    list_sm[[i]] <- tb %>% select(name, query_containment_ani) %>%
        mutate(genome_id = isolates$genome_id[i])
    # ANI estimated from the query containment in the match.
}
sm_genome <- bind_rows(list_sm) %>%
    select(genome_id, everything())
write_csv(sm_genome, paste0(folder_data, "genomics_analysis/taxonomy/sm_genome.csv"))
}

# 2. Aggregate the 16s blast results ----
# Read the reference
ref_16s_seq <- Biostrings::readDNAStringSet("~/bioinformatics/16s/refseq_16s.fasta")
ref_16s <- tibble(accession = str_remove(names(ref_16s_seq), "\\s.+"), scomment = str_remove(names(ref_16s_seq), "^[A-Z]+_[\\d]+.\\d\\s"))

# Read the blast results
list_b_16s <- rep(list(NA), nrow(isolates))
for (i in 1:nrow(isolates)) {
    list_b_16s[[i]] <- read_table(paste0(folder_data, "genomics/taxonomy/", isolates$genome_id[i],"/16s/blast_16s.txt"), col_names = names_blast) %>%
        mutate(genome_id = isolates$genome_id[i])
}

b_16s <- bind_rows(list_b_16s) %>%
    group_by(qseqid) %>%
    arrange(desc(bitscore)) %>%
    # Find the top hit
    slice(1) %>%
    mutate(accession = sseqid) %>%
    left_join(ref_16s) %>%
    select(genome_id, qseqid, scomment, everything()) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id)

write_csv(b_16s, paste0(folder_data, "genomics_analysis/taxonomy/b_16s.csv"))

# 3. Aggregate the genome blast results ----
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

b_genome <- bind_rows(list_b_genome) %>%
    #filter(pident > 90, bitscore > 10000, length > 10000) %>%
    filter(pident > 90, bitscore > 10000) %>%
    group_by(genome_id, qseqid) %>%
    arrange(desc(bitscore)) %>%
    # Find the top hit
    slice(1) %>%
    mutate(accession = sseqid) %>%
    left_join(ref_genome) %>%
    ungroup() %>%
    select(genome_id, species, strain, replicon, qseqid, scomment, sseqid, pident, length, bitscore) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id)

write_csv(ref_genome, paste0(folder_data, "genomics_analysis/taxonomy/ref_genome.csv"))
write_csv(b_genome, paste0(folder_data, "genomics_analysis/taxonomy/b_genome.csv"))

