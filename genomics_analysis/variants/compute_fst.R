#' This script computes Fst for each SNP

renv::load()
library(tidyverse)
library(apex) # for reading multiple fasta
library(poppr) # for processing fasta
library(mmod) # for computing Fst
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

#
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
list_sccg <- gpa$gene[which(apply(gpa[,-1], 1, sum) == 36)]
tb_sccg_fst <- tibble(gene = list_sccg, fst = NA)

for (i in 1:length(list_sccg)) {
    cat("\n",i , list_sccg[i])
    aln_file <- paste0(folder_data, "genomics/pangenome/isolates/aligned_gene_sequences/", list_sccg[i], ".aln.fas")
    fasta_data <- read.multiFASTA(aln_file)
    # Correct index names
    gid <- multidna2genind(fasta_data)
    indNames(gid) <- str_remove(indNames(gid), ";\\w+")
    # Assign populations
    isolates_pop <- tibble(genome_id = indNames(gid)) %>% left_join(isolates, by = join_by(genome_id))
    strata(gid) <- isolates_pop
    setPop(gid) <- ~site_group

    # Subset by pop
    cat("\televation")
    gid1 <- popsub(gid, c("high elevation", "low elevation"))
    vi1 <- diff_stats(gid1)$global
    cat("\turbanization")
    gid2 <- popsub(gid, c("suburban", "urban"))
    vi2 <- diff_stats(gid2)$global

    tb_sccg_fst$fst[i] <- list(tibble(metric = names(vi1), elev = vi1, urba = vi2))

}

sccg_fst <- tb_sccg_fst %>%
    #slice(1:2) %>%
    unnest(cols = fst)

write_csv(sccg_fst, paste0(folder_data, "genomics_analysis/variants/sccg_fst.csv"))

if (F) {
fasta_data <- read.multiFASTA(paste0(folder_data,
                                     "genomics/pangenome/isolates/aligned_gene_sequences/",
                                     "aarA", ".aln.fas"))
#plot(fasta_data)

# Correct index names
gid <- multidna2genind(fasta_data)
indNames(gid) <- str_remove(indNames(gid), ";\\w+")
# Assign populations
isolates_pop <- tibble(genome_id = indNames(gid)) %>% left_join(isolates)
strata(gid) <- isolates_pop
setPop(gid) <- ~site_group

# Subset by pop
gid_i <- popsub(gid, c("high elevation", "low elevation"))
xx <- diff_stats(gid_i)$global

a

gid_i <- popsub(gid, c("suburban", "urban"))
diff_stats(gid_i)
pairwise_Gst_Nei(gid_i, linearized = FALSE)

}

