#' This script aggregate the antiviral data

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv«"))
defense_finder <- read_csv(paste0(folder_data, "genomics_analysis/antiviral/defense_finder.csv"))
padloc <- read_csv(paste0(folder_data, "genomics_analysis/antiviral/padloc.csv"))

# meta data of the system database
padloc_meta <- read_delim(paste0(folder_data, "genomics_analysis/antiviral/padloc_meta.txt"), show_col_types = F) %>%
    rename(system_name = yaml.name)
    #mutate(system_namecleaned = str_replace_all(system, " ", "_") %>% tolower)


# Join isolate information
padloc <- padloc %>%
    distinct(genome_id, system.number, .keep_all = T) %>%
    left_join(select(isolates_tax, genome_id, contig_species)) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    rename(system_name = system) %>%
    #mutate(system_namecleaned = tolower(system)) %>%
    left_join(padloc_meta)

#
padloc %>%
    group_by(genome_id, contig_species) %>%
    count() %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = n), color = "black", fill = "grey90") +
    facet_grid2(.~contig_species, scales = "free_x", space = "free_x", strip = strip_vanilla(clip = "off")) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs()

p <- padloc %>%
    mutate(genome_id = rev(genome_id)) %>%
    group_by(contig_species, genome_id, group) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = group, y = genome_id, fill = n)) +
    facet_grid2(contig_species~., scales = "free_y", space = "free_y", switch = "y",
                strip = strip_themed(clip = "off", text_y = elem_list_text(angle = 0))) +
    coord_cartesian(clip = "off") +
    scale_fill_distiller(palette = 1, direction = 1) +
    theme_classic() +
    theme(
        panel.background = element_rect(color = NA, fill = "grey90"),
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = NA, fill = "white"),
        strip.placement = "outside"
    ) +
    guides() +
    labs(y = "")

ggsave(paste0(folder_data, "genomics_analysis/antiviral/", "01-padloc_count.png"), p, width = 6, height = 10)



