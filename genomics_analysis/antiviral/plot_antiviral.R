#' This script aggregate the antiviral data

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
defensefinder <- read_csv(paste0(folder_data, "genomics_analysis/antiviral/defensefinder.csv"))
padloc <- read_csv(paste0(folder_data, "genomics_analysis/antiviral/padloc.csv"))

# meta data of the system database
padloc_meta <- read_delim(paste0(folder_data, "genomics_analysis/antiviral/padloc_meta.txt"), show_col_types = F) %>%
    rename(system_name = yaml.name)

# Join isolate information
padloc <- padloc %>%
    distinct(genome_id, system.number, .keep_all = T) %>%
    left_join(select(isolates_tax, genome_id, contig_species)) %>%
    #mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    rename(system_name = system) %>%
    left_join(padloc_meta) %>%
    left_join(isolates) %>%
    mutate(genome_id = factor(genome_id))

defensefinder <- defensefinder %>%
    left_join(isolates) %>%
    left_join(select(isolates_tax, genome_id, contig_species)) %>%
    mutate(genome_id = factor(genome_id))

# 1. Padloc Count ----
plot_cols <- function (tb_cols) {
    tb_cols %>%
        ggplot() +
        geom_col(aes(x = genome_id, y = n), color = "black", fill = "grey90") +
        facet_grid2(contig_species~., scales = "free_y", space = "free_y", switch = "y",
                    strip = strip_themed(clip = "off", text_y = elem_list_text(angle = 0))) +
        coord_flip(clip = "off") +
        theme_classic() +
        theme(
            #panel.background = element_rect(color = NA, fill = "grey90"),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_rect(color = NA, fill = "white"),
            strip.placement = "outside"
        ) +
        guides() +
        labs(y = "")
}
plot_tiles <- function (tb_tiles, x) {
    tb_tiles %>%
        ggplot() +
        geom_tile(aes(x = {{x}}, y = genome_id, fill = n)) +
        facet_grid2(contig_species~., scales = "free_y", space = "free_y", switch = "y",
                    strip = strip_themed(clip = "off", text_y = elem_list_text(angle = 0))) +
        coord_cartesian(clip = "off") +
        scale_fill_distiller(palette = 1, direction = 1) +
        theme_classic() +
        theme(
            panel.background = element_rect(color = NA, fill = "grey90"),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            strip.placement = "outside"
        ) +
        guides() +
        labs(y = "")
}
p1 <- padloc %>%
    group_by(genome_id, contig_species) %>%
    count() %>%
    plot_cols

p2 <- padloc %>%
    group_by(contig_species, genome_id, system_name) %>%
    count() %>%
    plot_tiles(system_name) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

p <- plot_grid(p1, p2, nrow = 1, axis = "tb", align = "h", rel_widths = c(1,3))
ggsave(paste0(folder_data, "genomics_analysis/antiviral/", "01-padloc_count.png"), p, width = 20, height = 10)

# Padloc count by populatiosn in each gradient ----
plot_pops <- function (padloc_gra, x) {
    padloc_gra %>%
        group_by(population, contig_species, genome_id, {{x}}) %>%
        count() %>%
        ggplot() +
        geom_tile(aes(x = {{x}}, y = genome_id, fill = n)) +
        facet_grid2(population~., scales = "free_y", space = "free_y", switch = "y",
                    strip = strip_themed(clip = "off", text_y = elem_list_text(angle = 0))) +
        coord_cartesian(clip = "off") +
        scale_fill_distiller(palette = 1, direction = 1) +
        theme_classic() +
        theme(
            panel.background = element_rect(color = NA, fill = "grey90"),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            strip.placement = "outside"
        ) +
        guides() +
        labs(y = "")

}

p <- padloc %>%
    filter(contig_species == "S. medicae") %>%
    filter(gradient == "elevation") %>%
    plot_pops(group) +
    ggtitle("Elevation S. medicae")

ggsave(paste0(folder_data, "genomics_analysis/antiviral/", "elev_med/01-padloc_count.png"), p, width = 6, height = 5)

p <- padloc %>%
    filter(contig_species == "S. meliloti") %>%
    filter(gradient == "urbanization") %>%
    plot_pops(group) +
    ggtitle("Urbanization S. meliloti")
ggsave(paste0(folder_data, "genomics_analysis/antiviral/", "urbn_mel/01-padloc_count.png"), p, width = 6, height = 5)


# 2. Defense finder
p1 <- defensefinder %>%
    group_by(contig_species, genome_id) %>%
    count() %>%
    plot_cols
p2 <- defensefinder %>%
    group_by(contig_species, genome_id, type) %>%
    count() %>%
    plot_tiles(type) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

p <- plot_grid(p1, p2, nrow = 1, axis = "tb", align = "h", rel_widths = c(1,2))
ggsave(paste0(folder_data, "genomics_analysis/antiviral/", "02-defensefinder_count.png"), p, width = 15, height = 10)

##  defense finder count by populatiosn in each gradient ----
p <- defensefinder %>%
    filter(contig_species == "S. medicae") %>%
    filter(gradient == "elevation") %>%
    plot_pops(type) +
    ggtitle("Elevation S. medicae") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(paste0(folder_data, "genomics_analysis/antiviral/", "elev_med/02-defensefinder_count.png"), p, width = 6, height = 5)

p <- defensefinder %>%
    filter(contig_species == "S. meliloti") %>%
    filter(gradient == "urbanization") %>%
    plot_pops(type) +
    ggtitle("Urbanization S. meliloti") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
ggsave(paste0(folder_data, "genomics_analysis/antiviral/", "urbn_mel/02-defensefinder_count.png"), p, width = 6, height = 5)





