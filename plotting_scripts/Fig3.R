#' This script computes the effect size

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(lme4)
library(emmeans) # estimate marginal means
library(effectsize) # companion to Applied Regression
library(vcd) # for computing effect sizes of categorical response
library(boot)
source(here::here("metadata.R"))

set.seed(1)
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

compute_cohensd <- function (tb, response) {
    #' Compute Cohen's d
    formu <- paste0(response, " ~ site_group + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    emm.mod <- emmeans(mod, specs = "site_group")
    es <- eff_size(emm.mod, sigma = sigma(mod), edf = df.residual(mod))
    return(as_tibble(es))
}

list_treatments <- tibble(
    pop = c(rep("VA", 3), rep("PA", 3), rep("VA", 14), rep("PA", 7)),
    plant = c(rep("lupulina", 6), rep("sativa", 21)),
    nt = c(rep("without nitrogen", 6+7), rep("with nitrogen", 7), rep("without nitrogen", 7)),
    response = c(rep(c("nodule_number", "shoot_biomass_mg", "root_biomass_mg"), 2),
                 rep(c("nodule_number", "primary_root_nodule_number", "lateral_root_nodule_number", "shoot_height", "longest_petiole_length", "leaf_number", "leaf_color"), 2),
                 "nodule_number", "shoot_height", "leaf_number", "leaf_color", "primary_root_length", "lateral_root_number", "longest_lateral_root_length")
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
    mutate(
        pop = factor(pop, c("VA", "PA")),
        plant = factor(plant, c("lupulina", "sativa"))
    ) %>%
    mutate(response = str_remove(response, "mg") %>% str_replace_all("_", " ")) %>%
    arrange(pop, plant, response) %>%
    mutate(id = 1:n())

background_df <- tibble(pop = factor(rep(c("VA", "PA"), each = 2), c("VA", "PA")), plant = rep(c("sativa", "lupulina"), 2))


# VA population
p1 <- ess %>%
    filter(pop == "VA") %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_linerange(aes(x = response, ymin = lower_cl, ymax = upper_cl, shape = nt), color = "grey10", linewidth = 1, position = position_dodge2(width = .8)) +
    geom_point(aes(x = response, y = effect_size, shape = nt), size = 3, stroke = 1, fill = "white", position = position_dodge2(width = .8)) +
    scale_x_discrete(expand = c(0,.8), position = "top") +
    scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(plant~., space = "free_y", scales = "free_y", switch = "y") +
    coord_flip() +
    theme_bw() +
    theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.key = element_blank(),
        legend.key.spacing = unit(0, "mm"),
        legend.background = element_rect(color = NA, fill = "white"),
        plot.margin = margin(0,0,0,5, "mm"),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(color = "none", fill = "none") +
    labs(y = "standardized mean difference (Cohen's d)", title = "Elevation gradient")

# PA population
p2 <- ess %>%
    filter(pop == "PA") %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_linerange(aes(x = response, ymin = lower_cl, ymax = upper_cl, shape = nt), color = "grey10", linewidth = 1, position = position_dodge2(width = .8)) +
    geom_point(aes(x = response, y = effect_size, shape = nt), size = 3, stroke = 1, fill = "white", position = position_dodge2(width = .8)) +
    scale_x_discrete(expand = c(0,.8), position = "top") +
    scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(plant~., space = "free_y", scales = "free_y", switch = "y") +
    coord_flip() +
    theme_bw() +
    theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(0,0,0,5, "mm"),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(color = "none", fill = "none", shape = "none") +
    labs(y = "standardized mean difference (Cohen's d)", title = "Urbanization gradient")

p <- plot_grid(p1, p2, ncol = 1, scale = .95, align = "v", axis = "rl", labels = c("A", "B"), rel_heights = c(1.2,1)) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig3.png"), p, width = 5, height = 8)

