#' This script generates figure 2

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("analysis/00-metadata.R"))

# Read plant data
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
plants <- read_csv(paste0(folder_data, "temp/23-plants.csv"))
plants_long <- read_csv(paste0(folder_data, "temp/23-plants_long.csv"))
site_group_colors <- c(`high elevation` = "#0C6291", `low elevation` = "#BF4342", 
                 `suburban` = "#0cc45f", `urban` = "#a642bf", control = "grey")

# Read growth rate data
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'))

isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_longer(cols = -c(exp_id, temperature), names_to = "trait") %>%
    unite(trait, trait, temperature) %>%
    #pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD)) %>%
    left_join(isolates) 





# Panel A. Cartoons
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png")) + draw_text("placeholder")

# Panel B. Plot the 30C 
compute_trait_mean <- function (isolates_gc, tra = "r_30c", pop = "VA") {
igcl <- isolates_gc %>%
    filter(trait == tra) %>%
    filter(population == pop)
igcm <- igcl %>%
    group_by(population, site_group, trait) %>%
    summarize(mean = mean(value, na.rm = T), sem = sd(value, na.rm = T) / sqrt(n()))
    return(list(igcl = igcl, igcm = igcm))
}
plot_dots <- function (igcl, igcm) {
    set.seed(1)
    igcl %>%
        ggplot() +
        geom_rect(data = distinct(igcm, site_group), aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = igcm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = igcm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = unique(igcm$trait))
}

t2_1 <- compute_trait_mean(isolates_gc)
t2_2 <- compute_trait_mean(isolates_gc, pop = "PA")
igcl <- bind_rows(t2_1$igcl, t2_2$igcl) 
igcm <- bind_rows(t2_1$igcm, t2_2$igcm) 

# 
test_sign <- function (p) {
    if (p < 0.001) {
        x <- "***"
    } else if (p < 0.01) {
        x <- "**"
    } else if (p < 0.05) {
        x <- "*"
    } else {
        x <- "n.s."
    }
    return(x)
} 

# Does rhizobia sites have effect on shoot biomass?
isolates_test <- filter(isolates_gc, trait == "r_30c", population == "VA") %>% rename(r_30c = value)
mod <- lmer(r_30c ~ site_group + (1|site), data = isolates_test)
mod2_1 <- Anova(mod, type = 3) # no

# Does rhizobia sites have effect on shoot biomass?
isolates_test <- filter(isolates_gc, trait == "r_30c", population == "PA") %>% rename(r_30c = value)
mod <- lmer(r_30c ~ site_group + (1|site),, data = isolates_test)
mod2_2 <- Anova(mod, type = 3) # yes

bar_y <- 1.3
p2 <- igcl %>%
    ggplot() +
    geom_tile(data = igcm, aes(x = site_group, y = mean, fill = site_group), alpha = 0.3, height = Inf, width = 1) +
    geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", width = 0.1, height = 0) +
    geom_point(data = igcm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
    geom_errorbar(data = igcm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
    # Significance bars
    annotate("segment", x = 1, xend = 2, y = bar_y, yend = bar_y) +
    annotate("segment", x = 3, xend = 4, y = bar_y, yend = bar_y) +
    annotate("text", x = 1.5, y = bar_y, vjust = -1, label = test_sign(mod2_1[2,3])) + 
    annotate("text", x = 3.5, y = bar_y, vjust = -1, label = test_sign(mod2_2[2,3])) + 
    scale_fill_manual(values = site_group_colors) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(limits = c(0, bar_y*1.1)) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "r at 30C (1/hr)")


# Panel C. Plot the symbiosis traits comparing the two populations
compute_trait_mean2 <- function (plants_long, tra = "dry_weight_mg", pop = "VA") {
    pl <- plants_long %>%
        filter(trait == tra) %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots2 <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(data = distinct(plm, site_group), aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.2) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = "shoot biomass (mg)")
}

t3_1 <- compute_trait_mean2(plants_long,tra = "dry_weight_mg", pop = "VA")
t3_2 <- compute_trait_mean2(plants_long,tra = "dry_weight_mg", pop = "PA")
pl <- bind_rows(t3_1$pl, t3_2$pl) 
plm <- bind_rows(t3_1$plm, t3_2$plm) 

# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "VA", exp_id != "control")
mod <- lmer(dry_weight_mg ~ site_group + (1|plant) + (1|exp_id) + (1|waterblock), data = plants_test)
mod3_1 <- Anova(mod, type = 3) # no
mod <- lmer(nodule_number ~ site_group + (1|plant) + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # no
mod <- lmer(root_weight_mg ~ site_group + (1|plant) + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # no

# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "PA", exp_id != "control")
mod <- lmer(dry_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
mod3_2 <- Anova(mod, type = 3) # no
mod <- lmer(nodule_number ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes
mod <- lmer(root_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes

bar_y <- 55
p3 <- pl %>%
    ggplot() +
    geom_tile(data = plm, aes(x = site_group, y = mean, fill = site_group), alpha = 0.3, height = Inf, width = 1) +
    geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", width = 0.1, height = 0) +
    geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
    geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
    # Significance bars
    annotate("segment", x = 1, xend = 2, y = bar_y, yend = bar_y) +
    annotate("segment", x = 3, xend = 4, y = bar_y, yend = bar_y) +
    annotate("text", x = 1.5, y = bar_y, vjust = -1, label = test_sign(mod3_1[2,3])) + 
    annotate("text", x = 3.5, y = bar_y, vjust = -1, label = test_sign(mod3_2[2,3])) + 
    scale_fill_manual(values = site_group_colors) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(limits = c(0, bar_y*1.1)) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "shoot biomass (mg)")


# 
p_right <- plot_grid(p2, p3, nrow = 1, rel_heights = c(1, 1), scale = 0.95,
    align = "vh", axis = "tblr", labels = c("B", "C"))

p <- plot_grid(p1, p_right, rel_widths = c(1,1.5), labels = c("A", "")) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/Fig2.png"), p, width = 10, height = 4)


