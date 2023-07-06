#' Raw growth curve data

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L"))

# Subset only Ensifer
# subset_ensifer <- function(tb) {
#     tb %>%
#         left_join(select(isolates_RDP, strain, Genus)) %>%
#         drop_na()
# }
# gc <- gc %>% subset_ensifer()
# gc_summ <- gc_summ %>% subset_ensifer()
# gc.prm <- gc.prm %>% subset_ensifer()
# gc.prm.stat <- gc.prm.stat %>% subset_ensifer()


p <- gc_summ %>%
    mutate(strain = factor(strain, list_strains)) %>%
    left_join(rhizobia) %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs, color = rhizobia_site)) +
    geom_ribbon(aes(x = t, ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs, fill = rhizobia_site), alpha = 0.2) +
    facet_wrap(~strain, ncol = 5) +
    scale_x_continuous(breaks = seq(0, 48, 12)) +
    scale_color_manual(values = rhizobia_site_colors, labels = c("high", "low"), name = "site") +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), name = "site") +
    theme_classic() +
    theme(
        legend.position = c(0.9, 0.1),
        strip.background = element_rect(color = NA, fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(color = 1, fill = NA)
    ) +
    #guides(color = "none", fill = "none") +
    labs(x = "time (hrs)", y = expression(OD[600]))


ggsave(here::here("plots/FigS2.png"), p, width = 8, height = 6)
