#' This script plots the reaction norm of thermal adaptation
#' 1. Prepare the table
#' 2. Check model assumptions
#' 3. Run models
#' 4. Table
#' 5. Plot

library(tidyverse)
library(cowplot)
library(flextable)
library(ggh4x) # for nested facets
library(broom.mixed) # for tidying the model outputs
library(vegan)
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

# Prepare the data ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>% slice(1:32)
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))
gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(temperature, well, exp_id, r, lag, maxOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    drop_na(value)

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

# Stats ----
pairs_rn_perm <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_perm.csv"))
pairs_rn_posthoc <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_posthoc.csv"))

# Plot reaction norm ----
# Compute the mean
gtwlm <- gtwl %>%
    group_by(gradient, temperature, trait, population) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value))) %>%
    group_by(gradient, temperature, trait) %>%
    mutate(max_mean_value = max(mean_value, na.rm = T))

plot_rn <- function (gg, gra, tb_tidied) {
    #gra = "elevation"
    clean <- function (x) {
        x %>%
            mutate(trait = case_when(
            trait == "r" ~ "growth rate (1/hr)",
            trait == "lag" ~ "lag time (hr)",
            trait == "maxOD" ~ "yield [OD]"
        ))
    }
    # Mean value
    tb_mean <- filter(gtwlm, gradient == gra) %>%
        clean()

    # Stat
    tb_stat_perm <- pairs_rn_perm %>%
        rename(trait = trait_pre) %>%
        filter(gradient == gra) %>%
        mutate(ast = map_chr(p_value, turn_p_to_asteriks)) %>%
        filter(!str_detect(term, "Intercept")) %>%
        select(ii, gradient, trait, term, ast) %>%
        group_by(ii, gradient, trait) %>%
        pivot_wider(names_from = term, values_from = ast) %>%
        mutate(astlabs = paste0("population:temperature ", `population:temperature`, "\npopulation ", population, "\ntemperature ", temperature))
    tb_stat_posthoc <- pairs_rn_posthoc %>%
        rename(trait = trait_pre) %>%
        filter(gradient == gra) %>%
        mutate(ast = map_chr(p_value, turn_p_to_asteriks)) %>%
        select(ii, gradient, trait, temperature, ast) %>%
        left_join(select(tb_mean, gradient, temperature, trait, max_mean_value)) %>%
        filter(str_detect(ast, "\\*"))


    gg %>%
        filter(gradient == gra) %>%
        clean() %>%
        ggplot() +
        # Individual replicates
        geom_line(aes(x = temperature, y = value, group = well, color = population), alpha = 0.1) +
        # Mean value
        geom_ribbon(data = tb_mean, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = population, group = population), inherit.aes = FALSE, alpha = 0.2) +
        geom_point(data = tb_mean, aes(x = temperature, y = mean_value, color = population, group = population)) +
        geom_line(data = tb_mean, aes(x = temperature, y = mean_value, color = population, group = population)) +
        # Stats per panel
        geom_text(data = tb_stat_perm, aes(label = astlabs), x = 0.5, y = Inf, hjust = 0, vjust = 1.1, size = 2) +
        # stats per temperature
        geom_text(data = tb_stat_posthoc, aes(x = temperature, y = max_mean_value, label = ast), vjust = -.3, size = 5) +
        scale_color_manual(values = population_colors, name = "population") +
        scale_fill_manual(values = population_colors, name = "population") +
        scale_x_discrete(breaks = c("25c", "30c", "35c", "40c"), labels = c(25, 30, 35, 40)) +
        facet_wrap2(~trait, scales = "free_y", ncol = 1, strip.position = "left") +
        facetted_pos_scales(y = list(
                trait == "growth rate (1/hr)" ~ scale_y_continuous(breaks = seq(0, 1.5, .5), limits = c(0, 1.5)),
                trait == "lag time (hr)" ~ scale_y_continuous(breaks = seq(0, 60, 20), limits = c(0, 60)),
                trait == "yield [OD]" ~ scale_y_continuous(breaks = seq(0, .4, .2), limits = c(0, .45))
        )) +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey95", linewidth = .5),
            panel.border = element_rect(color = "black", fill = NA),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 10),
            strip.text.y = element_text(size = 10),
            strip.placement = "outside",
            axis.title.x = element_text(size = 10),
            axis.title.y = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.size = unit(5, "mm"),
            legend.box.margin = margin(0,0,0,0, "mm"),
            legend.margin = margin(0,0,0,0, "mm"),
            plot.background = element_blank(),
            plot.margin = margin(0,0,0,0, "mm")
        ) +
        guides(fill = guide_legend(override.aes = list(color = NA))) +
        labs(x = expression(paste("Temperature (", degree, "C)")))
}
p_rn1 <- plot_rn(gtwl, "elevation", tb_tidied)
p_rn2 <- plot_rn(gtwl, "urbanization", tb_tidied)

p <- plot_grid(
    p_rn1, p_rn2,
    nrow = 1, labels = c("A", "B"), scale = .95
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig4.png"), p, width = 6, height = 6)


# PERMANOVA ----
set.seed(1)
# Elevation
dat <- gtw %>%
    filter(temperature == "30c") %>%
    select(well, exp_id, r, lag, maxOD) %>%
    group_by(exp_id) %>%
    summarize(r = mean(r), lag = mean(lag), maxOD = mean(maxOD)) %>%
    left_join(select(isolates, population, exp_id, gradient)) %>%
    filter(gradient == "elevation")

pca_result <- select(dat, r, lag, maxOD) %>% prcomp(scale. = TRUE)
pcs <- as_tibble(pca_result$x) %>% mutate(population = dat$population)
dm <- vegdist(select(pcs, starts_with("PC")), method = "euclidean")
adonis2(dm ~ population, data = pcs, permutations = 1000)

# Urbnaization
dat <- gtw %>%
    filter(temperature == "30c") %>%
    select(well, exp_id, r, lag, maxOD) %>%
    group_by(exp_id) %>%
    summarize(r = mean(r), lag = mean(lag), maxOD = mean(maxOD)) %>%
    left_join(select(isolates, population, exp_id, gradient)) %>%
    filter(gradient == "urbanization")
pca_result <- select(dat, r, lag, maxOD) %>% prcomp(scale. = TRUE)
pcs <- as_tibble(pca_result$x) %>% mutate(population = dat$population)
dm <- vegdist(select(pcs, starts_with("PC")), method = "euclidean")
adonis2(dm ~ population, data = pcs, permutations = 1000)
