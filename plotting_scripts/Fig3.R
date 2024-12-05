#' This script plots the reaction norm of nitrogen treatments
#' 1. Prepare the table
#' 2. Plot
#' 3. permanova

library(tidyverse)
library(cowplot)
library(ggh4x) # for nested facets
library(vegan)
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

tidy_mod <- function (mod) {
    x <- tidy(mod)$estimate
    return(x)
}
boot_fun <- function(data, indices) {
    #' Define the bootstrapping function. This is generic depending on the exact model `mod` and the tidy function `tidy_mod`
    d <- data[indices, ]; eval(parse(text = st))
    #mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)
    return(tidy_mod(mod)) # Collect the fixed effect estimates
}
get_boot_cis <- function (boot_result) {
    #' Get the CIs of bootstrap value
    funn <- function(x) {
        tb <- as_tibble(x[4][[1]])
        colnames(tb) <- c("ci_lev", "ci_lower", "ci_upper")
        return(tb)
    }
    tibble(
        #term = tidy(mod)$term[1:4],
        t0 = boot_result$t0,
        ind = 1:length(boot_result$t0)
    ) %>%
        drop_na(t0) %>%
        mutate(
            conf = map(ind, ~boot.ci(boot_result, index = .x, type = "norm")),
            ci = map(conf, funn)
        ) %>%
        unnest(ci) %>%
        select(-conf) %>%
        select(ind, t0, everything())
}
do_stat <- function (dat, st) {
    d <- dat
    eval(parse(text = st))
    return(mod)
}

# 1. Prepare the data ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

plants_n <- plants %>%
    filter(population != "control", exp_plant == "sativa", gradient == "elevation") %>%
    select(
        -nodule_shape, -nodule_size, -nodule_color, -exp_labgroup,
        -primary_root_nodule_number, -lateral_root_nodule_number,
        -longest_petiole_length, -longest_lateral_root_length,
        -lateral_root_number, -primary_root_length
    ) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    ungroup()

plants_n_summ <- plants_n %>%
    group_by(exp_nitrogen, trait_type, trait_pre, population) %>%
    summarize(trait_mean = mean(value), trait_sem = sd(value)/sqrt(n())) %>%
    mutate(lower = trait_mean-qnorm(0.975)*trait_sem, upper = trait_mean+qnorm(0.975)*trait_sem)


# 2. Plot the reaction norm ----
p <- plants_n %>%
    group_by(gradient, population, exp_plant, exp_nitrogen, trait_type, trait_pre, value) %>%
    count() %>%
    ggplot(aes(x = exp_nitrogen, y = value)) +
    geom_point(aes(color = population, size = n), alpha = .5, shape = 16, position = position_dodge(width = .4)) +
    geom_linerange(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, ymin = lower, ymax = upper), linewidth = 1, position = position_dodge2(width = .4)) +
    geom_line(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, group = population), linewidth = 1, position = position_dodge(width = .4)) +
    geom_point(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean), size = 2, shape = 21, stroke = 1, fill = "white", position = position_dodge2(width = .4)) +
    scale_color_manual(values = population_colors) +
    scale_size_continuous(range = c(.5, 10)) +
    facet_nested(
        ~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none", nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1), solo_line = T,
        strip = strip_nested(bleed=T, clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("white", 3), rep("white", 7))))) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = NA)
    ) +
    guides(
        color = guide_legend(title = "population", override.aes = list(size = 5)),
        size = guide_legend(title = "sample size")
    ) +
    labs(x = "Nitrogen treatment", y = "")

ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 4)


# 3. PERMANOVA ----
set.seed(1)

dat <- plants %>%
    filter(population != "control", gradient == "elevation", exp_plant == "sativa") %>%
    drop_na(shoot_height, nodule_number, leaf_color, leaf_number) %>%
    select(gradient, population, site, exp_id, exp_nitrogen, shoot_height, nodule_number, leaf_color, leaf_number)
m <- select(dat, shoot_height, nodule_number, leaf_color, leaf_number)
adonis2(m ~ exp_nitrogen, data = dat, permutation = 1000)

