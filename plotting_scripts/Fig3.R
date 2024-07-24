#' This script

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
#library(ggh4x) # for nested facets
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
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
trait_names <- ess$response[c(1:3, 4,6,8,10,12,14:20)] %>%
    str_remove("_mg") %>%
    str_replace_all("_", " ")


# mountain population
p1 <- ess %>%
    filter(pop == "mountain") %>%
    mutate(id = c(1:3, rep(4:8, each = 2) + rep(c(-.1, .1), 5))) %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl), color = "grey10", linewidth = 1) +
    geom_point(aes(x = id, y = effect_size, shape = nt), size = 3, stroke = 1, fill = "white") +
    scale_x_continuous(breaks = 1:9, labels = trait_names[1:9], expand = c(0,.8)) +
    scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(.~plant, space = "free_x", scales = "free_x", switch = "y") +
    #coord_flip() +
    theme_bw() +
    theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        legend.key = element_blank(),
        legend.background = element_rect(color = "grey30", fill = "white"),
        #plot.margin = unit(c(0,0,0,0), "mm"),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(color = "none", fill = "none") +
    labs(y = "standardized mean difference (Cohen's d)")
# city population
p2 <- ess %>%
    filter(pop == "city") %>%
    mutate(id = c(1:3,4,5,7,8)) %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl), color = "grey10", linewidth = 1) +
    geom_point(aes(x = id, y = effect_size, shape = nt), size = 3, stroke = 1, fill = "white") +
    scale_x_continuous(breaks = 1:9, labels = trait_names[1:9], expand = c(0,.8)) +
    scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_fill_manual(values = plant_colors) +
    facet_grid(.~plant, space = "free_x", scales = "free_x", switch = "y") +
    #coord_flip() +
    theme_bw() +
    theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.1),
        legend.key = element_blank(),
        legend.background = element_rect(color = "grey30", fill = "white"),
        plot.margin = unit(c(10,10,10,0), "mm"),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(color = "none", fill = "none", shape = "none") +
    labs(y = "standardized mean difference (Cohen's d)")

p <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb", labels = c("A", "B")) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 5)

