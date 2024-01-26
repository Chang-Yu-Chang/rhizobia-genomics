#' This script plot the pangenome analysis

library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
library(vegan) # for dissimilarity matrix
library(ggtree) # for tree
library(ape) # for tree
source(here::here("analysis/00-metadata.R"))

# 0. read data ----
folder_roary <- paste0(folder_data, "temp/plasmidsaurus/summary/42-roary/roary5/")
# Read roary pangenome tables
gpa_meta <- read_csv(paste0(folder_roary, "gene_presence_absence.csv"), show_col_types = F) %>% clean_names()
gpa <- read_table(paste0(folder_roary, "gene_presence_absence.Rtab"), show_col_types = F)
gnew <- read_table(paste0(folder_roary, "number_of_new_genes.Rtab"), col_names = F, show_col_types = F)
gconserved <- read_table(paste0(folder_roary, "number_of_conserved_genes.Rtab"), col_names = F, show_col_types = F)
gunique <- read_table(paste0(folder_roary, "number_of_unique_genes.Rtab"), col_names = F, show_col_types = F)
gpan <- read_table(paste0(folder_roary, "number_of_genes_in_pan_genome.Rtab"), col_names = F, show_col_types = F)
# Read the mapping file
isolates_rhizo <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F)

# Remove rhizobium and keep only Ensifer
# gpa <- gpa %>%
#     select(-c(annotated_g1, annotated_g7, annotated_g12, annotated_g4, annotated_g18))

# pivot longer
gpa[gpa == 0] <- NA
gpa <- clean_names(gpa)
gpa_long <- gpa %>%
    pivot_longer(cols = -gene, values_drop_na = T, names_prefix = "annotated_", names_to = "genome_id") %>%
    mutate(genome_id = factor(genome_id, paste0("g", 1:19))) %>%
    clean_names()

# 0. some numbers ----
# Size of the Core Genome; gfs stands for gene frequency spectrum
gpa_glist <- gpa_long %>%
    group_by(gene) %>%
    count(name = "n_isolates") %>%
    mutate(n_isolates = factor(n_isolates, 1:19))

gpa_gfs <- gpa_glist %>%
    group_by(n_isolates) %>%
    count(name = "n_genes")
n_isolates <- nrow(gpa_gfs)
gpa_gfs$n_genes[gpa_gfs$n_isolates==n_isolates] # 2646 core genes when -i=75 and -cd=99 for 14 strains

# Pan Genome Size
length(unique(gpa_glist$gene)) # 17818

# Size of the Accessory Genome
gpa_gfs$n_genes

#  number of singleton
gpa_gfs$n_genes[1] # 7397


# 1. Number of gene shared in these isolates ----
#' This is the gene frequency spectrum G(k)
#' that correlates the number of orthologous genes groups (OGGs) containing genes from exactly k genomes
p <- gpa_gfs %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_col(aes(x = n_isolates, y = n_genes), color = "black", fill = "white") +
    #scale_x_continuous(breaks = 1:8) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(gpa_gfs$n_genes)/10^3*1.1), minor_breaks = 1:round(max(gpa_gfs$n_genes)/10^3,0)) +
    theme_bw() +
    theme() +
    guides() +
    labs(x = "gene shared by # of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-01-histogram.png"), plot = p, width = 5, height = 5)

# 2. heatmap for genes ----
p <- gpa_long %>%
    mutate(value = factor(value)) %>%
    ggplot() +
    geom_tile(aes(x = genome_id, y = gene, fill = value)) +
    scale_fill_manual(values = c(`1` = "maroon")) +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 10)
    ) +
    guides(fill = "none") +
    labs(x = "strain", y = "gene")
ggsave(paste0(folder_data, "temp/42-02-heatmap.png"), plot = p, width = 10, height = 20)


# 3. conserved gene vs. total genes in the pan genome ----
glist <- bind_rows(
    mutate(gconserved, key = "conserved genes", replicate = paste0("c", 1:10)),
    mutate(gpan, key = "total genes", replicate = paste0("t", 1:10))
) %>%
    pivot_longer(cols = -c(key, replicate), names_prefix = "X", names_to = "genome", values_to = "n_genes") %>%
    mutate(genome = as.numeric(genome), replicate = factor(replicate), key = factor(key, c("total genes", "conserved genes")))
p <- glist %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_line(aes(x = genome, y = n_genes, color = key, group = replicate), linewidth = .5) +
    scale_color_simpsons() +
    scale_x_continuous(minor_breaks = 1:n_isolates, limits = c(1,n_isolates)) +
    #scale_y_continuous(minor_breaks = 10, limits = c(1,20)) +
    theme_bw() +
    theme(
    ) +
    guides(color = guide_legend(title = "")) +
    labs(x = "# of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-03-conserved_vs_total_genes.png"), plot = p, width = 5, height = 3)

# numbers
range(glist$n_genes) # 2647 and 17818

# 4. new genes vs. unique genes ----
glist <- bind_rows(
    mutate(gunique, key = "unique genes", replicate = paste0("c", 1:10)),
    mutate(gnew, key = "new genes", replicate = paste0("t", 1:10))
) %>%
    pivot_longer(cols = -c(key, replicate), names_prefix = "X", names_to = "genome", values_to = "n_genes") %>%
    mutate(genome = as.numeric(genome), replicate = factor(replicate), key = factor(key, c("unique genes", "new genes")))
p <- glist %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_line(aes(x = genome, y = n_genes, color = key, group = replicate), linewidth = .5) +
    scale_color_simpsons() +
    scale_x_continuous(minor_breaks = 1:n_isolates, limits = c(1,n_isolates)) +
    #scale_y_continuous(minor_breaks = 10, limits = c(1,20)) +
    theme_bw() +
    theme(
    ) +
    guides(color = guide_legend(title = "")) +
    labs(x = "# of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-04-unique_vs_new_genes.png"), plot = p, width = 5, height = 3)

# Number
range(glist$n_genes) # 17 and 8418

# 5. all possible genes ----
glist <- bind_rows(
    mutate(gconserved, key = "conserved genes", replicate = paste0("c", 1:10)),
    mutate(gunique, key = "unique genes", replicate = paste0("u", 1:10)),
    mutate(gnew, key = "new genes", replicate = paste0("n", 1:10)),
    mutate(gpan, key = "total genes", replicate = paste0("t", 1:10)),
) %>%
    pivot_longer(cols = -c(key, replicate), names_prefix = "X", names_to = "genome", values_to = "n_genes") %>%
    mutate(genome = as.numeric(genome), replicate = factor(replicate), key = factor(key, c("total genes", "unique genes", "conserved genes", "new genes")))
p <- glist %>%
    mutate(n_genes = n_genes/1000) %>%
    ggplot() +
    geom_line(aes(x = genome, y = n_genes, color = key, group = replicate), linewidth = .5) +
    scale_color_simpsons() +
    scale_x_continuous(minor_breaks = 1:n_isolates, limits = c(1,n_isolates)) +
    theme_bw() +
    theme(
    ) +
    guides(color = guide_legend(title = "")) +
    labs(x = "# of genomes", y = expression("#" ~ of ~ genes ~ (10^3)))

ggsave(paste0(folder_data, "temp/42-05-all_genes.png"), plot = p, width = 5, height = 3)

# 6. build a tree using the gene presence and absence data ----
gpa_bin <- gpa
gpa_bin[is.na(gpa_bin)] <- 0
t_gpa <- gpa_bin %>%
    rename_with(~str_remove(., "^annotated_"), everything()) %>%
    select(-gene) %>%
    as.matrix() %>%
    t()
dim(t_gpa) # 14x17818
# Calculate Bray-Curtis dissimilarity
bin_dis <- vegdist(t_gpa, method = "bray")
clus <- hclust(bin_dis)

tree <- as.phylo(clus)
p <- tree %>%
    ggtree() +
    geom_tiplab()
ggsave(paste0(folder_data, "temp/42-06-tree_all_genes.png"), plot = p, width = 6, height = 3)

# 7. plot the key rhizobia genes ----
gpa_rhi <- gpa_long %>%
    filter(str_detect(gene, "nod|nif|fix")) %>%
    mutate(gene = str_sub(gene, 1, 4)) %>%
    distinct(gene, genome_id, value) %>%
    mutate(value = factor(value)) %>%
    mutate(genome_id = factor(genome_id, paste0("g", 20:1))) %>%
    mutate(gene_group = str_sub(gene, 1, 3))

p1 <- gpa_rhi %>%
    # add a pseudo block for plotting absence data
    #bind_rows(tibble(gene = "fixX", genome_id = "g2", value = "0", gene_group = "fix")) %>%
    mutate(genome_id = factor(genome_id, paste0("g", 20:1))) %>%
    mutate(value = factor(value, c("1", "0"))) %>%
    ggplot() +
    geom_tile(aes(x = gene, y = genome_id, fill = value), color = "black", width = 0.8, height = 0.8, linewidth = .8) +
    facet_grid(~gene_group, scales = "free_x", space = "free") +
    scale_fill_manual(values = c(`1` = "grey20", `0` = "white"), labels = c(`1` = "presence", `0` = "absence")) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    # scale_x_discrete(position = "top", expand = c(0,0)) +
    # scale_y_discrete(position = "left", expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 10),
        strip.text = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
        #legend.position = "bottom"
    ) +
    guides(fill = "none") +
    #guides(fill = guide_legend(title = NULL, override.aes = list(color = "black", linewidth = .5))) +
    labs(x = "gene", y = "genome")

# check the growth curve trait, scaled
if (T) {
    # Call the isolates data. THIS IS A DIRTY CODE FOR CONVIENIENCE
    gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
    gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
    gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
    gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)
    isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
        rename(strain = ExpID) %>%
        filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L"))

    # Subset only Ensifer
    subset_ensifer <- function(tb) {
        tb %>%
            left_join(select(isolates_RDP, strain, Genus)) %>%
            drop_na()
    }

    gc <- gc %>% subset_ensifer()
    gc_summ <- gc_summ %>% subset_ensifer()
    gc.prm <- gc.prm %>% subset_ensifer()
    gc.prm.stat <- gc.prm.stat %>% subset_ensifer()
    isolates_rhizo <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F)

    isolates <- isolates_rhizo %>%
        select(strain = exp_id, genus, genome_id) %>%
        left_join(gc.prm.stat) %>%
        mutate(genome_id = factor(genome_id, paste0("g", 1:19))) %>%
        filter(genus == "Ensifer")

    write_csv(isolates, paste0(folder_data, "temp/42-isolates.csv"))
    # Clean this for anvio input data layer
    isolates_anvio <- isolates %>%
        select(genome_id, r, lag, maxOD) %>%
        mutate(genome_id = str_replace(genome_id, "g", "Chang_Q5C_")) %>%
        rename(sample = genome_id) %>%
        mutate(r_category = ifelse(r > median(r), "high_r", "low_r"))
    write_delim(isolates_anvio, paste0(folder_data, "temp/42-isolates_anvio.txt"), delim = "\t")

    "CLEAN THIS ISOLATES TRAIT DATA IN A BETTER FORM"

}
isolates_long <- isolates %>%
    select(genome_id, r, lag, maxOD) %>%
    # scale
    mutate_at(c("r", "lag", "maxOD"), ~(scale(.) %>% as.vector)) %>%
    pivot_longer(cols = -genome_id, names_to = "trait")

p2 <- isolates_long %>%
    mutate(genome_id = factor(genome_id, paste0("g", 20:1))) %>%
    mutate(trait = factor(trait, c("r", "lag", "maxOD"))) %>%
    ggplot() +
    geom_tile(aes(x = trait, y = genome_id, fill = value), color = "black", width = 0.8, height = 0.8, linewidth = .8) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "maroon", breaks = c(1, 0, -1), na.value = "white") +
    scale_color_manual(values = c(g19="white", g2="black", g4="black", g8="black", g10="black", g13="black", g15="black")) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    # geom_tile(aes(x = trait, y = genome_id, fill = value)) +
    # geom_vline(xintercept = c(1:10)-0.5, color = "black") +
    # geom_hline(yintercept = c(1:20)-0.5, color = "black") +
    #facet_grid(~., scales = "free_x", space = "free") +
    # scale_fill_gradient2(low = "steelblue", mid = "white", high = "maroon", breaks = c(1, 0, -1)) +
    # scale_x_discrete(position = "top", expand = c(0,0)) +
    # scale_y_discrete(position = "left", expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none"
    ) +
    guides(fill = guide_legend(title = NULL, override.aes = list(color = "black", linewidth = .5))) +
    labs(x = "growth trait")

# Plot the plant phenotypes
# gc_labels <- gc.prm.stat %>%
#     mutate(strain_label = factor(1:n())) %>%
#     select(strain, strain_label)
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))

isolates_plant <- treatments_M %>%
    group_by(strain) %>%
    summarize(mean_nodule_number = mean(nodule_number, na.rm = T),
              mean_dry_weight_mg = mean(dry_weight_mg, na.rm = T),
              mean_root_weight_mg = mean(root_weight_mg, na.rm = T)) %>%
    left_join(isolates) %>%
    select(genome_id, starts_with("mean")) %>%
    filter(!is.na(genome_id))

isolates_long3 <- isolates_plant %>%
    # scale
    mutate_at(c("mean_nodule_number", "mean_dry_weight_mg", "mean_root_weight_mg"), ~(scale(.) %>% as.vector)) %>%
    pivot_longer(cols = -genome_id, names_to = "trait")

p3 <- isolates_long3 %>%
    bind_rows(tibble(genome_id = "g19", trait = "mean_dry_weight_mg", value = NA)) %>%
    mutate(genome_id = factor(genome_id, paste0("g", c(19,17:15, 13, 11:8, 6:2)))) %>%
    # shorten trait names
    mutate(trait = str_remove(trait, "mean_") %>% str_remove("_weight_mg")) %>%
    mutate(trait = str_replace(trait, "nodule_number", "#nodule") %>% str_replace("dry", "shoot")) %>%
    ggplot() +
    geom_tile(aes(x = trait, y = genome_id, fill = value, color = genome_id),  width = 0.8, height = 0.8, linewidth = .8) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "maroon", breaks = c(1, 0, -1), na.value = "white") +
    scale_color_manual(values = c(g19="white", g2="black", g4="black", g8="black", g10="black", g13="black", g15="black")) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left", drop = F) +
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        #panel.background = element_rect(fill = "grey", color = NA),
        legend.position = "right"
    ) +
    guides(color ="none" ,fill = guide_legend(title = "standardized\ntrait value", override.aes = list(color = "black", linewidth = .5))) +
    labs(x = "symbiosis trait")

p <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 0.23, 0.38))

ggsave(paste0(folder_data, "temp/42-07-rhizobia_genes.png"), plot = p, width = 12, height = 6)




















