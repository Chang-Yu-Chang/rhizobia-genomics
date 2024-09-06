#' This script

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
# library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for permanova
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>% slice(1:32)
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

set.seed(1)
detect_sig <- function (pv) {
    if (pv > 0.05) {
        return("n.s.")
    } else if (pv > 0.01) {
        return("*")
    } else if (pv > 0.001) {
        return("**")
    } else if (pv < 0.001) {
        return("***")
    }
}

isolates <- isolates %>%
    mutate(population = case_when(
        population == "VA" ~ "elevation",
        population == "PA" ~ "urbanization"
    ))

gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T))

# PCA at 30C ----
gts1 <- gtwl %>%
    filter(population == "elevation", temperature == "30c") %>%
    pivot_wider(names_from = trait, values_from = value) %>%
    select(site_group, genome_id, r, lag, maxOD)
gts2 <- gtwl %>%
    filter(population == "urbanization", temperature == "30c") %>%
    pivot_wider(names_from = trait, values_from = value) %>%
    select(site_group, genome_id, r, lag, maxOD)

pca_results <- list(elevation = prcomp(gts1[,-c(1,2)], scale. = TRUE),
                    urbanization = prcomp(gts2[,-c(1,2)], scale. = TRUE))
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

pcs <- bind_rows(
    as_tibble(pca_results$elevation$x) %>% mutate(site_group = gts1$site_group, genome_id = gts1$genome_id),
    as_tibble(pca_results$urbanization$x) %>% mutate(site_group = gts2$site_group, genome_id = gts2$genome_id)
) %>%
    left_join(distinct(isolates, population, site_group)) %>%
    mutate(population = factor(population, c("elevation", "urbanization")))

plot_pca <- function (pcs, pop) {
    pcsi <- pcs %>% filter(population == pop)
    dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
    mod <- with(pcsi, adonis2(dm ~ site_group, data = pcsi, permutations = 10000, strata = genome_id))
    pcsi %>%
        ggplot() +
        geom_point(aes(x = PC1, y = PC2, color = site_group), shape = 21, stroke = 1, size = 2) +
        stat_ellipse(aes(x = PC1, y = PC2, fill = site_group), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
        geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
        geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
        scale_color_manual(values = site_group_colors) +
        scale_fill_manual(values = site_group_colors, name = "population") +
        scale_x_continuous(breaks = seq(-8,8,2)) +
        scale_y_continuous(breaks = seq(-8,8,2)) +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey10", fill = NA),
            plot.background = element_blank(),
            legend.background = element_blank()
        ) +
        guides(color = "none") +
        labs(x = paste0("PC1 (", get_pcvar(pca_results[pop][[1]])[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results[pop][[1]])[2], "%)"))

}
p_pca1 <- plot_pca(pcs, "elevation")
p_pca2 <- plot_pca(pcs, "urbanization")

p_pca <- plot_grid(p_pca1, p_pca2, nrow = 1, axis = "lr", align = "v", labels = c("B", "C"))

# Stats
pcsi <- pcs %>% filter(population == "elevation")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- with(pcsi, adonis2(dm ~ site_group, data = pcsi, permutations = 10000, strata = genome_id))

pcsi <- pcs %>% filter(population == "urbanization")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- with(pcsi, adonis2(dm ~ site_group, data = pcsi, permutations = 10000, strata = genome_id))



# Reaction norm ----
gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    mutate(trait = case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "maxOD" ~ "yield [OD]"
    ))

# Compute the mean
gtwlm <- gtwl %>%
    group_by(population, temperature, trait, site_group) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value)))

plot_rn <- function (gg, pop) {
    tt <- filter(gtwlm, population == pop)

    gg %>%
        filter(population == pop) %>%
        ggplot() +
        geom_line(aes(x = temperature, y = value, group = well, color = site_group), alpha = 0.1) +
        # Mean value
        geom_ribbon(data = tt, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = site_group, group = site_group), inherit.aes = FALSE, alpha = 0.2) +
        geom_point(data = tt, aes(x = temperature, y = mean_value, color = site_group, group = site_group)) +
        geom_line(data = tt, aes(x = temperature, y = mean_value, color = site_group, group = site_group)) +
        scale_color_manual(values = site_group_colors, name = "population") +
        scale_fill_manual(values = site_group_colors, name = "population") +
        facet_wrap(~trait, scales = "free_y", nrow = 1, strip.position = "left") +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 10),
            strip.text.y = element_text(size = 10),
            strip.placement = "outside",
            axis.title.x = element_text(size = 10),
            axis.title.y = element_blank(),
            legend.position = "none",
            plot.background = element_blank()
        ) +
        guides() +
        labs()

}

p_rn1 <- plot_rn(gtwl, "elevation")
p_rn2 <- plot_rn(gtwl, "urbanization")
p_rn <- plot_grid(p_rn1, NULL, p_rn2, nrow = 1, axis = "lr", align = "v", rel_widths = c(1,.1,1), labels = c("D", "", "E"))


# Combine figures ----

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig2.png"), scale = 1) +
    draw_plot(p_pca, x = .3, y = .48, width = .7, height = .45) +
    # draw_plot(p_pca1, x = .3, y = .48 , width = .35, height = .45) +
    # draw_plot(p_pca2, x = .65, y = .48, width = .35, height = .45) +
    draw_plot(p_rn, x = .05, y = .05, width = .9, height = .4) +
    # draw_plot(p_rn1, x = .01, y = .05, width = .45, height = .4) +
    # draw_plot(p_rn2, x = .5, y = .05, width = .45, height = .4) +
    #draw_plot(p_sativa, x = .35, y = -0.01, width = .6, height = .32) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig2.png"), p, width = 12, height = 6)


if (F) {

    # r_30c ~ elevation w/ N
    tb <- gts %>%
        left_join(isolates) %>%
        mutate(population = factor(population, c("VA", "PA"))) %>%
        filter(population == "VA")
    mod <- lmer(r_30c ~ site_group + (1|site), data = tb)
    Anova(mod, type = 3)
    vmax <- max(tb$r_30c)*.95

    p1 <- tb %>%
        ggplot() +
        geom_boxplot(aes(x = site_group, y = r_30c, fill = site_group)) +
        geom_jitter(aes(x = site_group, y = r_30c, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
        annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
        annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
        scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
        scale_color_manual(values = alpha(site_group_colors, 1)) +
        #scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1),
            axis.title.y = element_text(size = 10),
            panel.border = element_rect(color = "grey10", fill = NA),
            plot.background = element_blank()
        ) +
        guides(fill = "none", color = "none") +
        labs(y = "growth rate at 30C (1/hr)")


    tb <- gts %>%
        left_join(isolates) %>%
        mutate(population = factor(population, c("VA", "PA"))) %>%
        filter(population == "PA")
    mod <- lmer(r_30c ~ site_group + (1|site), data = tb)
    Anova(mod, type = 3)
    vmax <- max(tb$r_30c)*.95

    p2 <- tb %>%
        ggplot() +
        geom_boxplot(aes(x = site_group, y = r_30c, fill = site_group)) +
        geom_jitter(aes(x = site_group, y = r_30c, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
        annotate("segment", x = 1.05, xend = 1.95, y = vmax, yend = vmax) +
        annotate("text", x = 1.5, y = vmax*1.05, label = detect_sig(Anova(mod)[3])) +
        scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
        scale_color_manual(values = alpha(site_group_colors, 1)) +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1),
            axis.title.y = element_text(size = 10),
            panel.border = element_rect(color = "grey10", fill = NA),
            plot.background = element_blank()
        ) +
        guides(fill = "none", color = "none") +
        labs(y = "growth rate at 30C (1/hr)")
}
