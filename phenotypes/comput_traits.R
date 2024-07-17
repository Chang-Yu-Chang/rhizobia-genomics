#' This script computes the rhizobia traits data

renv::load()
library(tidyverse)
library(cowplot)
library(ggh4x) # for nested facets
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for computing jaccard and adonis2
source("metadata.R")

plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
gts <- read_csv(paste0(folder_phenotypes, 'growth/gts.csv'))

isolates_gc <- gts %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_longer(cols = -c(exp_id, temperature), names_to = "trait") %>%
    unite(trait, trait, temperature) %>%
    left_join(isolates)

background_df <- tibble(population = c("VA", "VA", "PA", "PA"), site_group = names(site_group_colors)[1:4])

# Plot the boxplot
p <- isolates_gc %>%
    filter(str_detect(trait, "30c")) %>%
    #filter(str_detect(trait, "r_30c")) %>%
    #filter(trait %in% c("r_30c", "")) %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_jitter(aes(x = site_group, y = value), alpha = 0.5, shape = 16, size = 3, color = "black", width = 0.1, height = 0) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_fill_manual(values = site_group_colors) +
    facet_nested(trait~population + site_group, scales = "free", switch = "y") +
    #facet_grid(trait~population, scales = "free", switch = "y") +
    theme_minimal() +
    theme(
        strip.placement = "outside",
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing.x = unit(c(0,10,0), "mm"),
        panel.spacing.y = unit(5, "mm"),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(fill = "none") +
    labs(x = "")

ggsave(paste0(folder_phenotypes, "traits/01-growth.png"), p, width = 5, height = 8)







