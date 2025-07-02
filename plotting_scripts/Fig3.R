#' GCV

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
library(lme4)
library(car)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
tt <- read_gpas()
tb_cog <- tibble(genome_id = factor(iso$genome_id, iso$genome_id)) %>%
    left_join(select(iso, genome_id, population, contig_species)) %>%
    mutate(cog_classify = map(genome_id, ~read_tsv(paste0(folder_data, "genomics/cog/", .x, "/cog_classify.tsv")))) %>%
    unnest(cog_classify) %>%
    clean_names()


# Panel A. gcv  ----
p1 <- tbtr$tr[[2]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = 0, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    #scale_x_continuous(limits = c(0, 100)) +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 4, y = 25) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()


p1_1 <- isolates %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    select(genome_id, population, contig_species) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    ggplot() +
    geom_tile(aes(x = population, y = genome_id, fill = contig_species), color = "black", linewidth = .5) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors, breaks = c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.text = element_text(face = "italic", size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.background = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
    labs()


# Panel B. symb genes ----
tb <- tt$gd %>%
    #filter(!str_detect(gene, "group")) %>%
    filter(str_detect(gene, "nod|nif|fix|noe|exo|nol|fdx")) %>%
    mutate(ge = str_sub(gene, 1, 5) %>% str_remove("\\d$|_$")) %>%
    mutate(g = str_sub(ge, 1, 3)) %>%
    select(g, ge, gene, genome_id) %>%
    left_join(select(iso, genome_id, contig_species))

p2 <- tb %>%
    mutate(
        genome_id = factor(genome_id, rev(get_taxa_name(p1))),
        contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens", "control"))
    ) %>%
    group_by(g, contig_species, genome_id, ge) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = ge, y = genome_id, fill = n)) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "grey80", high = "grey20") +
    facet_grid2(contig_species ~ g, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        #axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        legend.position = "right",
        strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.background.x = element_rect(color = NA, fill = "gray90")
    ) +
    guides() +
    labs()

## Stat: do s meliloti and s medicate differ in their genes?
# mod <- tb %>%
#     filter(contig_species %in% c("S. medicae", "S. meliloti")) %>%
#     select(ge, genome_id, contig_species) %>%
#     group_by(ge, genome_id, contig_species) %>%
#     count() %>%
#     lmer(n ~ ge + contig_species + (1|contig_species:genome_id), data = .)
#Anova(mod, type = 3)
#emmeans(mod, ~contig_species)

# Panel C. cog ----
p3 <- tb_cog %>%
    mutate(
        genome_id = factor(genome_id, rev(get_taxa_name(p1))),
        contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens", "control"))
    ) %>%
    group_by(cog_letter, cog_id, contig_species, genome_id) %>%
    filter(cog_letter %in% LETTERS[20:26]) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = cog_id, y = genome_id), fill = "maroon", color = NA) +
    coord_cartesian(clip = "off") +
    facet_grid2(~cog_letter, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        panel.spacing.x = unit(0, "cm")
    ) +
    guides() +
    labs()



# ----
p <- plot_grid(
    p1, p1_1 + guides(fill = "none"),
    p2, p3, nrow = 1,
    align = "h", axis = "tb", rel_widths = c(1,.12,1.5, 1.5),
    labels = c("A", "", "B", "C")
) +
    draw_text("GCV", x = .22, y = .9, size = 10, hjust = 0) +
    draw_plot(get_legend(p1_1), x = -.41, y = .29) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 6)

