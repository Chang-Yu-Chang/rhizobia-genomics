#' COG

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

iso <- read_csv(paste0(folder_data, "output/iso.csv"))
symbiosis_genes <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/symbiosis_genes.csv"))
tt <- read_gpas()

#
tb <- tibble(genome_id = factor(iso$genome_id, iso$genome_id)) %>%
    left_join(select(iso, genome_id, population, contig_species)) %>%
    mutate(cog_classify = map(genome_id, ~read_tsv(paste0(folder_data, "genomics/cog/", .x, "/cog_classify.tsv")))) %>%
    unnest(cog_classify) %>%
    clean_names()

# Find the core and acceer gernes of each sp
s1 <- iso$genome_id[iso$contig_species == "S. meliloti"]
temp <- tb %>%
    filter(genome_id %in% s1) %>%
    distinct(cog_id, genome_id) %>%
    group_by(cog_id) %>%
    count()
list_core1 <- temp$cog_id[temp$n == length(s1)]
list_acce1 <- temp$cog_id[temp$n < length(s1)]

#
s2 <- iso$genome_id[iso$contig_species == "S. medicae"]
temp <- tb %>%
    filter(genome_id %in% s2) %>%
    distinct(cog_id, genome_id) %>%
    group_by(cog_id) %>%
    count()
list_core2 <- temp$cog_id[temp$n == length(s2)]
list_acce2 <- temp$cog_id[temp$n < length(s2)]


#
x1 <- tb %>%
    filter(contig_species == "S. meliloti") %>%
    select(population, contig_species, genome_id, cog_letter, cog_id) %>%
    filter(cog_id %in% list_core1) %>%
    group_by(population, contig_species, genome_id, cog_letter) %>%
    distinct(genome_id, cog_id) %>%
    count()

x2 <- tb %>%
    filter(contig_species == "S. meliloti") %>%
    select(population, contig_species, genome_id, cog_letter, cog_id) %>%
    filter(cog_id %in% list_acce1) %>%
    group_by(population, contig_species, genome_id, cog_letter) %>%
    count()

x3 <- tb %>%
    filter(contig_species == "S. medicae") %>%
    select(population, contig_species, genome_id, cog_letter, cog_id) %>%
    filter(cog_id %in% list_core2) %>%
    group_by(population, contig_species, genome_id, cog_letter) %>%
    count()

x4 <- tb %>%
    filter(contig_species == "S. medicae") %>%
    select(population, contig_species, genome_id, cog_letter, cog_id) %>%
    filter(cog_id %in% list_acce2) %>%
    group_by(population, contig_species, genome_id, cog_letter) %>%
    count()

p <- bind_rows(
    mutate(x1, gene_group = "core"),
    mutate(x2, gene_group = "acce"),
    mutate(x3, gene_group = "core"),
    mutate(x4, gene_group = "acce")
) %>%
    filter(gene_group == "core") %>%
    ggplot() +
    #geom_boxplot(aes(x = contig_species, y = n, color = population), position = position_dodge(width = .8)) +
    geom_point(aes(x = contig_species, y = n, color = population), shape = 21, position = position_jitterdodge(dodge.width = .5, jitter.height = 0, jitter.width = .1)) +
    facet_wrap(~cog_letter+gene_group, scales = "free_y", nrow = 4) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(here::here("plots/FigS6.png"), p, width = 12, height = 8)

