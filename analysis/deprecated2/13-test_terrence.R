#' This script cleans Terrence's data

library(tidyverse)
library(cowplot)
library(janitor)
library(ggsci)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))

treatments_ttb <- readxl::read_xlsx(paste0(folder_data, "raw/rhizobia/Terrence/wood Lab rotation -fall 22 .xlsx"))
treatments_cyc <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    clean_names()

# Clean Terrence's data ----
treatments_ttb <- treatments_ttb %>%
    rename(id = plant, waterblock = block,
           rhizobia = rhizobia_strain, rhizobia_site = location,
           dry_weight_mg = shoot_weight, nodule_number = nodules, root_weight_mg = root_weight) %>%
    # Correct unit
    mutate(dry_weight_mg = dry_weight_mg*10^3, root_weight_mg = root_weight_mg*10^3) %>%
    # Plant unique id
    mutate(id = id + nrow(treatments_cyc)) %>%
    mutate(rhizobia_site = tolower(rhizobia_site)) %>%
    arrange(id) %>%
    clean_names()


# 0.2 merge the data ----
treatments <- treatments_cyc %>% bind_rows(treatments_ttb) %>%
    replace_na(list(rhizobia = "control", rhizobia_site = "control")) %>%
    mutate(rhizobia = str_replace(rhizobia, "control_\\d", "control")) %>%
    mutate(rhizobia_site = case_when(
        rhizobia_site == "H" ~ "high-elevation",
        rhizobia_site == "L" ~ "low-elevation",
        T ~ rhizobia_site
    )) %>%
    mutate(rhizobia_site = factor(rhizobia_site, c("high-elevation", "low-elevation", "urban", "suburban", "control")))

nrow(treatments) # 251 plants
treatments %>% filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>% nrow # 231 plants with non-zero shoot mass

# 1. plot all plant traits ----
treatments_long <- treatments %>%
    select(id, rhizobia_site, rhizobia, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>%
    pivot_longer(cols = c(-id, -rhizobia, -rhizobia_site), names_to = "trait")
length(unique(treatments_long$id)) # 231 plants
length(unique(treatments_long$rhizobia)) # 14 rhizobia + 1 control
treatments_long %>%
    arrange(rhizobia_site) %>%
    distinct(id, rhizobia_site, rhizobia) %>%
    group_by(rhizobia_site, rhizobia) %>%
    count()


p <- treatments_long %>%
    ggplot() +
    geom_boxplot(aes(x = rhizobia, y = value, color = rhizobia), outlier.shape = NA) +
    geom_jitter(aes(x = rhizobia, y = value, color = rhizobia), width = 0.2, height = 0, shape = 21) +
    facet_grid(trait~rhizobia_site, scales = "free", space = "free_x") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_rect(fill = NA, color = NA)
    ) +
    guides(color = "none") +
    labs()

ggsave(paste0(folder_data, "temp/13-01-plant_traits.png"), p, width = 10, height = 10)

# 2. Plot the nodule vs. plant biomass ----
treatments_wide <- treatments %>%
    select(id, rhizobia_site, rhizobia, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0)

p <- treatments_wide %>%
    filter(rhizobia_site != "control") %>%
    ggplot() +
    geom_point(aes(x = nodule_number, y = dry_weight_mg, color = rhizobia_site), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
    scale_color_manual(values = brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)]) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    coord_equal() +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey95")
    ) +
    guides() +
    labs(x = "# of nodules", y = "shoot biomass (mg)")

ggsave(paste0(folder_data, "temp/13-02-nodule_vs_shoot.png"), p, width = 5, height = 4)

# 3. facet plot
p <- treatments_wide %>%
    filter(rhizobia_site != "control") %>%
    ggplot() +
    geom_point(aes(x = nodule_number, y = dry_weight_mg, color = rhizobia_site), shape = 21, size = 2, stroke = 1)+
    scale_color_manual(values = brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)]) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    facet_wrap(~rhizobia_site, nrow = 2) +
    coord_equal() +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey95"),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = "# of nodules", y = "shoot biomass (mg)")

ggsave(paste0(folder_data, "temp/13-03-nodule_vs_shoot_facet.png"), p, width = 8, height = 6)

# Stat for slope in regression
lm_x_y <- function (dat, x_var, y_var) {
    formula <- as.formula(paste(y_var, "~", x_var))
    lm(formula, data = dat) %>%
        broom::tidy() %>%
        filter(term == x_var) %>%
        select(term, estimate, p.value)
}

treatments_wide %>%
    filter(rhizobia_site != "control") %>%
    nest(data = c(id, rhizobia, dry_weight_mg, nodule_number, root_weight_mg)) %>%
    mutate(lm = map(data, ~lm_x_y(.x, "nodule_number", "dry_weight_mg"))) %>%
    select(rhizobia_site, lm) %>%
    unnest(lm)










