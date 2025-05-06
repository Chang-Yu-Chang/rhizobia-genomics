#'

library(tidyverse)
library(cowplot)
library(waffle)
source(here::here("metadata.R"))

# Read data
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    mutate(contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")))
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>% filter(population != "mid elevation") %>%
    # Use sites where the isolates were from
    filter(site %in% isolates$site)
dml <- read_csv(paste0(folder_phenotypes, "sites/dml.csv"))
tb_month <- read_csv(paste0(folder_phenotypes, "sites/tb_month.csv"))

# Panel A. site temperature ----
tb_tmax <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    group_by(population) %>%
    summarize(mean_tmax = mean(tmax_deg_c))

p1 <- dml %>%
    filter(yday >= 182 & yday <= 273) %>%
    ggplot() +
    geom_histogram(aes(x = tmax_deg_c), position = "identity", alpha = .5) +
    geom_vline(data = tb_tmax, aes(xintercept = mean_tmax)) +
    facet_grid(~population) +
    coord_flip(clip = "off") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 10),
        strip.background = element_blank()
    ) +
    guides() +
    labs(x = expression("Daily maximum temperature "(degree*C)), y = "Num. of days in Jul-Sep")

# Panel B. strains ----
p2 <- iso %>%
    group_by(population, contig_species) %>%
    count() %>%
    ggplot(aes(fill = contig_species, values = n)) +
    geom_waffle(n_rows = 2, flip = T, color = "white", radius = unit(1, "mm")) +
    scale_fill_manual(values = species_colors) +
    scale_y_continuous(breaks = c(5, 10)) +
    facet_grid(~population) +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(10, "mm"),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0, 12, 0, 15), unit = "mm")
    ) +
    guides(fill = guide_legend(nrow = 2)) +
    labs()


p <- plot_grid(p1, p2, ncol = 1) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig1.png"), p, width = 3, height = 6)

