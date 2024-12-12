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

# 1. Prepare the data ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
nitrogen_rn_perm <- read_csv(paste0(folder_phenotypes, "nitrogen_rn/nitrogen_rn_perm.csv"))

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
    pivot_longer(cols = -c(1:12), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    ungroup()

plants_n_summ <- plants_n %>%
    group_by(exp_nitrogen, trait_type, trait_pre, population) %>%
    summarize(trait_mean = mean(value), trait_sem = sd(value)/sqrt(n())) %>%
    mutate(lower = trait_mean-qnorm(0.975)*trait_sem, upper = trait_mean+qnorm(0.975)*trait_sem) %>%
    mutate(trait_type = factor(trait_type, unique(traits$trait_type)), trait_pre = factor(trait_pre, traits$trait_pre))


# 2. Plot the reaction norm ----
# Stat
tb_stat_perm <- nitrogen_rn_perm %>%
    mutate(ast = map_chr(p_value, turn_p_to_asteriks)) %>%
    filter(!str_detect(term, "Intercept")) %>%
    select(ii, trait_type, trait_pre, term, ast) %>%
    group_by(ii, trait_type, trait_pre) %>%
    pivot_wider(names_from = term, values_from = ast) %>%
    mutate(trait_type = factor(trait_type, unique(traits$trait_type)), trait_pre = factor(trait_pre, traits$trait_pre)) %>%
    mutate(astlabs = paste0("population:nitrogen ", `population:exp_nitrogen`, "\npopulation ", population, "\nnitrogen ", exp_nitrogen))

p <- plants_n %>%
    group_by(gradient, population, exp_plant, exp_nitrogen, trait_type, trait_pre, value) %>%
    mutate(trait_type = factor(trait_type, unique(traits$trait_type)), trait_pre = factor(trait_pre, traits$trait_pre)) %>%
    count() %>%
    ggplot(aes(x = exp_nitrogen, y = value)) +
    geom_point(aes(color = population, size = n), alpha = .5, shape = 16, position = position_dodge(width = .4)) +
    geom_linerange(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, ymin = lower, ymax = upper), linewidth = 1, position = position_dodge2(width = .4)) +
    geom_line(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, group = population), linewidth = 1, position = position_dodge(width = .4)) +
    geom_point(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean), size = 2, shape = 21, stroke = 1, fill = "white", position = position_dodge2(width = .4)) +
    # Stats per panel
    geom_text(data = tb_stat_perm, aes(label = astlabs), x = 0.5, y = Inf, hjust = 0, vjust = 1.1, size = 2) +
    scale_color_manual(values = population_colors) +
    scale_size_continuous(range = c(.5, 10)) +
    scale_y_continuous(expand = c(0.1, 2)) +
    facet_nested(
        ~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none", nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1), solo_line = T,
        strip = strip_nested(bleed=T , clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("white", 3), rep("white", 7))))) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = NA)
    ) +
    guides(
        color = guide_legend(title = "population", override.aes = list(size = 5)),
        size = guide_legend(title = "num. of plants")
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

