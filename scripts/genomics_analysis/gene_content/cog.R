#'

library(tidyverse)
library(janitor)
library(lme4)
library(car)
library(emmeans)
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

bind_rows(
    mutate(x1, gene_group = "core"),
    mutate(x2, gene_group = "acce"),
    mutate(x3, gene_group = "core"),
    mutate(x4, gene_group = "acce")
) %>%
    filter(gene_group == "core") %>%
    ggplot() +
    geom_jitter(aes(x = contig_species, y = n, color = population), width = .1, shape = 21, height = 0) +
    facet_grid(~cog_letter+gene_group, scales = "free_y", nrow = 1) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()




if (F) {

#  stat
tbc <- tb %>%
    select(population, contig_species, genome_id, cog_letter) %>%
    group_by(population, contig_species, genome_id, cog_letter) %>%
    filter(contig_species %in% c("S. medicae", "S. meliloti")) %>%
    count()


mod <- lmer(n ~ cog_letter + population + contig_species + (1|genome_id), data = tbc)
Anova(mod, type = 3)

emm <- emmeans(mod, ~contig_species)
pai <- pairs(emm) %>% as_tibble()

pai %>%
    separate(contrast, into = c("a", "b"), sep = " - ") %>%
    mutate(cog1 = str_sub(a, str_count(a), str_count(a))) %>%
    mutate(cog2 = str_sub(b, str_count(b), str_count(b))) %>%
    filter(cog1 == cog2) %>%
    filter(str_detect(a, "medicae"), str_detect(b, "meliloti"))

}
#
tb %>%
    select(population, contig_species, genome_id, cog_letter, cog_id) %>%
    group_by(population, contig_species, genome_id, cog_letter, cog_id) %>%
    count()


tb %>%
    group_by(population, contig_species, genome_id, cog_letter) %>%
    count() %>%
    #mutate(n = log(n)) %>%
    filter(cog_letter %in% LETTERS[21:26]) %>%
    ggplot() +
    #geom_tile(aes(x = cog_letter, y = genome_id, fill = n)) +
    geom_col(aes(x = genome_id, y = n)) +
    facet_grid(cog_letter~contig_species, scales = "free", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()


tb %>%
    filter(cog_letter == "X") %>%
    view


# bar plot
tb %>%
    group_by(contig_species, genome_id, cog_letter) %>%
    #filter(cog_id == "COG0071") %>%
    #filter(str_detect(cog_name, " heat ")) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = n, fill = cog_letter), color = 1) +
    facet_grid(~contig_species, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()

# col by letter

tb %>%
    group_by(contig_species, genome_id, cog_letter) %>%
    #filter(cog_id == "COG0071") %>%
    #filter(str_detect(cog_name, " heat ")) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = n), color = 1) +
    facet_grid(cog_letter~contig_species, scales = "free", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()



# Tile by cog id
sort(unique(tb$cog_id))
tb %>%
    #filter(genome_id %in% c("g2", "g3")) %>%
    #filter(cog_id %in% paste0("COG", sprintf("%04d", 1:10))) %>%
    group_by(contig_species, genome_id, cog_letter,cog_id) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = cog_id, y = genome_id, fill = n)) +
    facet_grid(contig_species~cog_letter, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    guides() +
    labs()

# Tile by cog letter
sort(unique(tb$cog_id))
tb %>%
    group_by(contig_species, genome_id, cog_letter) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = cog_letter, y = genome_id, fill = log(n))) +
    geom_text(aes(x = cog_letter, y = genome_id, label = n)) +
    facet_grid(contig_species~., scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.ticks.x = element_blank()
    ) +
    guides() +
    labs()


tb %>%
    group_by(population, contig_species, genome_id, cog_letter) %>%
    filter(contig_species %in% c("S. medicae", "S. meliloti")) %>%
    count() %>%
    ggplot() +
    #geom_boxplot(aes(x = cog_letter, y = n)) +
    geom_jitter(aes(x = population, y = n, color = contig_species, shape = population), width = .2, stroke = .8, alpha = .8) +
    scale_shape_manual(values = c(VA = 21, PA = 22)) +
    scale_color_manual(values = species_colors) +
    #geom_text(aes(x = cog_letter, y = genome_id, label = n)) +
    facet_grid(~cog_letter, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.ticks.x = element_blank()
    ) +
    guides() +
    labs()

# tree
tbw <- tb %>%
    group_by(genome_id, cog_id) %>%
    count() %>%
    pivot_wider(names_from = cog_id, values_from = n, values_fill = 0)

# zz <- tbw[,-1]
# rownames(zz) <- tbw$genome_id
# tr <- dist(zz) %>% hclust(method = "ward.D2") %>% ape::as.phylo()
