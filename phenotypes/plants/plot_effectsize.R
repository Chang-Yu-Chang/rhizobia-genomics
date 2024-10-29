#' Plot effect size

library(tidyverse)
source(here::here("metadata.R"))

# 1. Cohen's d
cohensds <- read_csv(paste0(folder_data, "phenotypes/effectsize/cohensds.csv"))
plot_cohensds <- function (cohensds, plant) {
    cohensds %>%
        filter(exp_plant == plant) %>%
        mutate(exp_nitrogen = factor(exp_nitrogen, c("without nitrogen", "with nitrogen"))) %>%
        ggplot() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_linerange(aes(x = trait, ymin = lower_cl, ymax = upper_cl), color = "grey10", linewidth = 1, position = position_dodge2(width = .5)) +
        geom_point(aes(x = trait, y = effect_size, shape = exp_nitrogen), size = 3, stroke = 1, fill = "white", position = position_dodge2(width = .5, reverse = T)) +
        scale_x_discrete(expand = c(0,.8), position = "top") +
        scale_y_continuous(limits = c(-3, 3), expand = c(0,.1), breaks = -3:3) +
        scale_shape_manual(values = c(`with nitrogen` = 16, `without nitrogen` = 21), labels = c("N+", "N-")) +
        scale_fill_manual(values = plant_colors) +
        facet_grid(gradient~exp_nitrogen, scales = "free_y", space = "free_y", switch = "y") +
        coord_flip(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            plot.background = element_rect(color = NA, fill = "white")
        ) +
        guides(color = "none", fill = "none") +
        labs(y = "standardized mean difference (Cohen's d)", title = paste0(plant, " experiment"))
}
p <- plot_cohensds(cohensds, "sativa")
ggsave(paste0(folder_phenotypes, "effectsize/01-cohensd_sativa.png"), p, width = 8, height = 5)
p <- plot_cohensds(cohensds, "lupulina")
ggsave(paste0(folder_phenotypes, "effectsize/01-cohensd_lupulina.png"), p, width = 5, height = 5)






if (F) {


    # 2. Plot Hedge's g ----
    ess <- list_treatments %>%
        unnest(hedgesg) %>%
        clean_names() %>%
        mutate(id = 1:n()) %>%
        mutate(nt = case_when(
            nt == "without nitrogen" ~ "w/o N",
            nt == "with nitrogen" ~ "w/ N"
        )) %>%
        mutate(pop = case_when(
            pop == "PA" ~ "urbanization",
            pop == "VA" ~ "elevation"
        )) %>%
        mutate(study_name = paste0(" ", response))

    p <- ess %>%
        ggplot() +
        geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl), color = "grey10", linewidth = 1) +
        geom_point(aes(x = id, y = effect_size, shape = nt), size = 3, stroke = 1) +
        scale_x_reverse(breaks = 1:20, labels = ess$study_name, expand = c(0,.8), position = "top") +
        scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
        scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
        scale_fill_manual(values = plant_colors) +
        coord_flip() +
        facet_nested(pop + plant ~., scale = "free_y", space = "free_y", switch = "y", strip = strip) +
        theme_minimal() +
        theme(
            axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
            panel.grid.minor.y = element_blank(),
            panel.spacing.y = unit(0, "mm"),
            strip.text = element_text(size = 10),
            plot.background = element_rect(color = NA, fill = "white"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.background = element_rect(color = "black", fill = "white")
        ) +
        guides(color = "none", fill = "none") +
        labs(y = "standardized mean difference (Hedge's g)")

    ggsave(paste0(folder_phenotypes, "effectsize/02-hedgesg.png"), p, width = 8, height = 6)
# 3. Plot eta squared partial ----
ess <- list_treatments %>%
    unnest(partialetasquared) %>%
    clean_names() %>%
    mutate(id = 1:n()) %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(study_name = paste(nt, response, sep = ", "))

p <- ess %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = id, y = eta2_partial, color = plant), size = 3) +
    geom_segment(aes(x = id, xend = id, y = ci_low, yend = ci_high, color = plant), linewidth = 1) +
    scale_x_reverse(breaks = 1:13, labels = ess$study_name, expand = c(0,.7)) +
    scale_color_aaas() +
    coord_flip() +
    facet_grid(plant ~., scale = "free_y", space = "free_y") +
    theme_light() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 15)
    ) +
    guides(color = "none") +
    labs(y = "partial eta squared")
ggsave(paste0(folder_phenotypes, "effectsize/03-partialetasquared.png"), p, width = 8, height = 6)




}

