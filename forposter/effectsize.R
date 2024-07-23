#' This script computes the effect sizes of plant experiments

renv::load()
library(tidyverse)
library(janitor)
library(ggsci)
library(ggh4x) # for nested facets
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(emmeans) # estimate marginal means
library(effectsize) # companion to Applied Regression
library(vcd) # for computing effect sizes of categorical response
library(boot)
source(here::here("metadata.R"))

#
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
lupulinas <- read_csv(paste0(folder_phenotypes, "plants/lupulinas.csv"))
sativas <- read_csv(paste0(folder_phenotypes, "plants/sativas.csv"))

# Remove control and nonsymbiontic strains
plants <- plants %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")
lupulinas <- lupulinas %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")
sativas <- sativas %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    select(genome_id = genome, species = sm_species, exp_id) %>%
    bind_rows(tibble(genome_id = c("g_src1", "g_bg1"), species = NA, exp_id = c("src-1", "bg-1")))

# Cohen's d
compute_cohensd <- function (tb, response) {
    formu <- paste0(response, " ~ site_group + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    emm.mod <- emmeans(mod, specs = "site_group")
    es <- eff_size(emm.mod, sigma = sigma(mod), edf = df.residual(mod))
    return(as_tibble(es))
}


#
list_treatments <- tibble(
    pop = c(rep("VA", 3), rep("PA", 3), rep("VA", 10), rep("PA", 4)),
    plant = c(rep("lupulina", 6), rep("sativa", 14)),
    nt = c(rep("without nitrogen", 6+5), rep("with nitrogen", 5), rep("without nitrogen", 4)),
    response = c(rep(c("nodule_number", "shoot_biomass_mg", "root_biomass_mg"), 2),
                 rep(c("nodule_number", "shoot_height", "longest_petiole_length", "leaf_number", "leaf_color"), 2),
                 "nodule_number", "shoot_height", "leaf_number", "leaf_color")
)

list_treatments <- list_treatments %>%
    rowwise() %>%
    mutate(tb = list(get(paste0(plant, "s")) %>% filter(exp_nitrogen == nt, population == pop))) %>%
    mutate(cohensd = list(compute_cohensd(tb, response)))

# Plot cohen's d
ess <- list_treatments %>%
    unnest(cohensd) %>%
    clean_names() %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(pop = case_when(
        pop == "PA" ~ "city",
        pop == "VA" ~ "mountain"
    )) %>%
    mutate(
        study_name = paste0(" ", response),
        pop = factor(pop, c("mountain", "city")),
        plant = factor(plant, c("lupulina", "sativa"))
    ) %>%
    arrange(pop, plant, response) %>%
    mutate(id = 1:n())

background_df <- tibble(pop = factor(rep(c("mountain", "city"), each = 2), c("mountain", "city")), plant = rep(c("sativa", "lupulina"), 2))
strip <- strip_nested(background_y = list(
    NULL, NULL,
    element_rect(fill = alpha(plant_colors[2], 0.2), color = NA),
    element_rect(fill = alpha(plant_colors[1], 0.2), color = NA),
    element_rect(fill = alpha(plant_colors[2], 0.2), color = NA),
    element_rect(fill = alpha(plant_colors[1], 0.2), color = NA)
))


ess$study_name[c(1:3, 4,6,8,10,12,14:20)]

trait_names <- ess$response[c(1:3, 4,6,8,10,12,14:20)] %>%
    str_remove("_mg") %>%
    str_replace_all("_", " ")

p <- ess %>%
    mutate(id = c(1:3, rep(4:8, each = 2) + rep(c(-.1, .1), 5) , 9:15)) %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl), color = "grey10", linewidth = 1) +
    geom_point(aes(x = id, y = effect_size, shape = nt), size = 3, stroke = 1) +
    scale_x_continuous(breaks = 1:15, labels = trait_names, expand = c(0,.8), position = "top") +
    scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_fill_manual(values = plant_colors) +
    coord_flip() +
    facet_nested(pop + plant~., scale = "free_y", space = "free_y", strip = strip, switch = "y") +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.spacing.y = unit(c(0,5,0), "mm"),
        strip.text = element_text(size = 10),
        plot.background = element_rect(color = NA, fill = "white"),
        legend.position = c(0.85, 0.6),
        legend.title = element_blank(),
        legend.background = element_rect(color = "grey30", fill = "grey90")
    ) +
    guides(color = "none", fill = "none") +
    labs(y = "standardized mean difference (Cohen's d)")

ggsave(here::here("forposter/effectsize.pdf"), p, width = 8, height = 5)


