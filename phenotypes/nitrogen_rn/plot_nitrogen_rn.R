#' This script plots the reaction norm of nitrogen treatments

library(tidyverse)
library(flextable)
library(cowplot)
library(ggh4x) # for nested facets
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
nitrogen_rn_perm <- read_csv(paste0(folder_phenotypes, "nitrogen_rn/nitrogen_rn_perm.csv"))

# Plot the reaction norm ----
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

p <- plants_n %>%
    group_by(gradient, population, exp_plant, exp_nitrogen, trait_type, trait_pre, value) %>%
    count() %>%
    ggplot(aes(x = exp_nitrogen, y = value)) +
    geom_point(aes(color = population, size = n), alpha = .5, shape = 16, position = position_dodge(width = .3)) +
    geom_linerange(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, ymin = lower, ymax = upper), linewidth = 1, position = position_dodge2(width = .3)) +
    geom_line(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, group = population), linewidth = 1, position = position_dodge(width = .3)) +
    geom_point(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean), size = 2, shape = 21, stroke = 1, fill = "white", position = position_dodge2(width = .3)) +
    scale_color_manual(values = population_colors) +
    scale_size_continuous(range = c(.5, 10)) +
    facet_nested(
        ~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none", nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1), solo_line = T,
        strip = strip_nested(bleed=T, clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("white", 3), rep("white", 7))))) +
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

ggsave(paste0(folder_phenotypes, "nitrogen_rn/nitrogen_rn.png"), p, width = 10, height = 4)

# Make the table ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: ")) %>%
        str_remove(fixed(": value")) %>%
        str_replace("exp_id", "strain") %>%
        str_replace("exp_labgroup", "labgroup")
}
ft2 <- nitrogen_rn_perm  %>%
    select(ii, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Estimate = statistic, P = siglab) %>%
    # Clean the table
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre)
    ) %>%
    select(-ii) %>%
    arrange(Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Type", "Trait", "Model")) %>%
    valign(j = c("Type", "Trait", "Model"), valign = "center") %>%
    align(j = c("Type", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Model", space = 1.5) %>%
    width(j = "Model", width = 4) %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft2, path = paste0(folder_phenotypes, "nitrogen_rn/nitrogen_rn_perm.png"), res = 300)
