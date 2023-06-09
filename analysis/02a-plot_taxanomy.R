#' This scripts plots the isolated strain taxonomy

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))


isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    filter(!is.na(Family))

# 1. isolate taxonomy by family and by genus ----
p1 <- isolates_RDP %>%
    ggplot() +
    geom_bar(aes(y = Family, fill = Family), color = 1) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(isolates_RDP)), vjust = 2, hjust = 2) +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = grey(0.8), linetype = 1))

p2 <- isolates_RDP %>%
    ggplot() +
    geom_bar(aes(y = Genus, fill = Family), color = 1) +
    annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(isolates_RDP)), vjust = 2, hjust = 2) +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = grey(0.8), linetype = 1))

p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")

ggsave(paste0(folder_data, "temp/02a-01-isolate_taxonomy.png"), plot = p, width = 8, height = 7)

# 2. Isolate ID distribution ----
isolates_ID <- isolates_ID %>%
    mutate(Site = ExpID %>% str_sub(1, 3) %>% str_replace("-", "") %>% str_replace("M", "") %>%
               str_replace("fp1|fp2", "fp") %>% str_replace("gp1", "gp")) %>%
    mutate(Owner = ifelse(str_sub(Site, 1, 1) %in% c("L", "H"), "CYC", "TPP"))

site_id <- unique(isolates_ID$Site)
family_id <- c("Rhizobiaceae", "Pseudomonadaceae", "Enterobacteriaceae", "Others")


p <- isolates_RDP %>%
    mutate(Family = as.character(Family)) %>%
    mutate(Family = ifelse(Family %in% c("Rhizobiaceae", "Pseudomonadaceae", "Enterobacteriaceae"), Family, "Others")) %>%
    left_join(isolates_ID) %>%
    mutate(Site = factor(Site, site_id)) %>%
    mutate(Family = factor(Family, family_id)) %>%
    group_by(Owner, Site, Family) %>%
    summarize(Count = n()) %>%
    ggplot() +
    geom_col(aes(x = Site, y = Count, fill = Family), color = 1) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(breaks = 1:10) +
    theme_classic()
ggsave(paste0(folder_data, "temp/02a-02-family_count.png"), plot = p, width = 6, height = 3)



