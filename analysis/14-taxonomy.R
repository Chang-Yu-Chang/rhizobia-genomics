#' This script assigns the taxonomy 

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
isolates <- isolates %>% drop_na(exp_id)
names_blast <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# 1. Aggregate the 16s blast results
# Read the reference 
ref_16s_seq <- Biostrings::readDNAStringSet("~/bioinformatics/16s/refseq_16s.fasta")
ref_16s <- tibble(accession = str_remove(names(ref_16s_seq), "\\s.+"), 
    scomment = str_remove(names(ref_16s_seq), "^[A-Z]+_[\\d]+.\\d\\s")) 

# Read the blast results 
list_b_16s <- rep(list(NA), nrow(isolates))
for (i in 1:nrow(isolates)) {
    list_b_16s[[i]] <- read_table(paste0(folder_data, "genomics/taxonomy/", isolates$genome_id[i],"/16s/blast_16s.txt"), col_names = names_blast) %>%
        mutate(genome_id = isolates$genome_id[i])
}

b_16s <- bind_rows(list_b_16s)

b_16s <- b_16s %>%
    group_by(qseqid) %>%
    arrange(desc(bitscore)) %>%
    slice(1) %>%
    mutate(accession = sseqid) %>%
    left_join(ref_16s) %>%
    select(genome_id, qseqid, scomment, everything()) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id)

write_csv(b_16s, paste0(folder_data, "temp/14-b_16s.csv"))

# 2. Aggregate the genome blast results
ref_genome_seq <- Biostrings::readDNAStringSet(paste0(folder_data, "genomics/blast_db/genomes.fasta"))
ref_genome <- tibble(accession = str_remove(names(ref_genome_seq), "\\s.+"), 
    scomment = str_remove(names(ref_genome_seq), "^[A-Z]+_[A-Z\\d]+.\\d\\s")) 
ref_genome <- ref_genome %>%
    rowwise() %>%
    mutate(species = str_split(scomment, " ")[[1]] %>% `[`(2)) %>%
    mutate(strain = str_remove(scomment, "Sinorhizobium \\w+ |Ensifer \\w+ ") %>% str_remove("strain ") %>% str_remove(",")) %>%
    mutate(strain = str_split(strain, " ")[[1]] %>% `[`(1)) %>%
    mutate(replicon = case_when(
        str_detect(scomment, "chromosome") ~ "chromosome",
        str_detect(scomment, "ctg") ~ str_extract(scomment, "ctg\\d+"),
        #str_detect(scomment, "ctg") ~ str_replace(scomment, ".+ctg", "ctg", ) %>% str_replace("ctg\\d+ .+", "ctg\\d+"),
        str_detect(scomment, "plasmid") ~ str_extract(scomment, "plasmid [A-Z|0-9|a-z|_]+")
    )) 
# Manually correct the comments
# NC_003047.1 is 1021 chromsome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006965.1/
# NC_009636.1 is WSM419 chromosome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000017145.1/
# NC_012587.1 is NGR234 chromosome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000018545.1/
ref_genome$replicon[ref_genome$accession == "NC_003047.1"] <- "chromosome"
ref_genome$replicon[ref_genome$accession == "NC_009636.1"] <- "chromosome"
ref_genome$replicon[ref_genome$accession == "NC_012587.1"] <- "chromosome"

# Read the blast results
list_b_genome <- rep(list(NA), nrow(isolates))
for (i in 1:nrow(isolates)) {
    list_b_genome[[i]] <- read_table(paste0(folder_data, "genomics/taxonomy/", isolates$genome_id[i],"/blast_genome/blast_genome.txt"), col_names = names_blast) %>%
        mutate(genome_id = isolates$genome_id[i])
}

b_genome <- bind_rows(list_b_genome)

b_genome <- b_genome %>%
    filter(pident > 90, bitscore > 10000) %>%
    group_by(genome_id, qseqid) %>%
    arrange(desc(bitscore)) %>%
    #slice(1) %>%
    mutate(accession = sseqid) %>%
    left_join(ref_genome) %>%
    select(genome_id, qseqid, scomment, everything()) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id)

write_csv(ref_genome, paste0(folder_data, "temp/14-ref_genome.csv"))
write_csv(b_genome, paste0(folder_data, "temp/14-b_genome.csv"))

# Out put for each contig, the best alignment
isolates_contigs <- read_csv(paste0(folder_data, 'temp/12-contigs.csv'))

contigs <- b_genome %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    mutate(contig_id = paste0(genome_id, "_", qseqid)) %>%
    select(genome_id, contig_id, species, strain, replicon, bitscore) %>%
    group_by(genome_id, contig_id) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(mutate(isolates_contigs, contig_id = paste0(genome_id, "_", contig_id))) %>%
    mutate(replicon = case_when(
        str_detect(replicon, "psymA") ~ str_replace(replicon, "psymA", "pSymA"),
        str_detect(replicon, "psymB") ~ str_replace(replicon, "psymB", "pSymB"),
        str_detect(replicon, "pA") ~ str_replace(replicon, "pA", "pSymA"),
        str_detect(replicon, "pB") ~ str_replace(replicon, "pB", "pSymB"),
        str_detect(replicon, "pSMED01") ~ str_replace(replicon, "pSMED01", "pSymA like"),
        str_detect(replicon, "pSMED02") ~ str_replace(replicon, "pSMED02", "pSymB like"),
        str_detect(replicon, "pSMED03") ~ str_replace(replicon, "pSMED03", "accessory"),
        str_detect(replicon, "pWSM1115_1") ~ str_replace(replicon, "pWSM1115_1", "pSymA like"),
        str_detect(replicon, "pWSM1115_2") ~ str_replace(replicon, "pWSM1115_2", "pSymB like"),
        str_detect(replicon, "pWSM1115_3") ~ str_replace(replicon, "pWSM1115_3", "accessory"),
        str_detect(replicon, "pSU277_1") ~ str_replace(replicon, "pSU277_1", "pSymA like"),
        str_detect(replicon, "pSU277_2") ~ str_replace(replicon, "pSU277_2", "pSymB like"),
        str_detect(replicon, "pSU277_3") ~ str_replace(replicon, "pSU277_3", "accessory"),
        str_detect(replicon, "pSU277_4") ~ str_replace(replicon, "pSU277_4", "accessory"),
        T ~ replicon
    )) %>%
    filter(!str_detect(replicon, "ctg"))
unique(contigs$replicon)
#view(ref_genome)
write_csv(contigs, paste0(folder_data, "temp/14-contigs.csv"))

# 3. Assign isolates to taxonomy
contigs <- read_csv(paste0(folder_data, "temp/12-contigs.csv"))
isolates_contigs <- b_genome %>%
    rename(contig_id = qseqid) %>%
    left_join(contigs) %>%
    filter(contig_length > 1000000) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id, desc(contig_length)) %>%
    ungroup() 

write_csv(isolates_contigs, paste0(folder_data, "temp/14-isolates_contigs.csv"))

# Addison's strains
isolates_abm <- isolates_contigs %>%
    filter(genome_id %in% paste0("g", 38:43)) %>%
    left_join(select(isolates, genome_name, genome_id)) %>%
    select(genome_name, genome_id, contig_id, contig_length, scomment, everything()) %>%
    group_by(genome_id, contig_id) %>%
    slice(1:10)
write_csv(isolates_abm, paste0(folder_data, "temp/14-isolates_abm.csv"))
