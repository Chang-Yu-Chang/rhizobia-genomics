#' Growth traits for all strains

library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("analysis/00-metadata.R"))

gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L"))

gc_labels <- gc.prm.stat %>%
    mutate(strain_label = factor(1:n())) %>%
    select(strain, strain_label)

# By strain
plot_strains_trait <- function (gc.prm, ytrait, ylab = "") {
    gc.prm %>%
        left_join(gc_labels) %>%
        #select(-t.r, - startOD) %>%
        #pivot_longer(cols = c(r, lag, maxOD), names_to = "trait") %>%
        ggplot() +
        geom_rect(data = tibble(strain_site_group = c("H", "L")), aes(fill = strain_site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_point(aes(x = strain_label, y = {{ytrait}}, group = strain_label, color = strain_label), shape = 21, size = 2, stroke = 1, fill = NA,
                   position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
        scale_color_manual(values = rep("black", 100)) +
        scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
        facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high elevation", L="low elevation"))) +
        theme_classic() +
        theme(
            panel.grid.major.x = element_line(color = "grey80"),
            panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
            strip.background = element_rect(color = NA, fill = NA),
            strip.text = element_text(size = 10, color = "black"),
            strip.placement = "outside",
            axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 10, color = "black"),
            legend.position = "none",
            plot.margin = unit(c(0,5,0,0), "mm")
        ) +
        guides(color = "none") +
        labs(x = "rhizobia strain", y = ylab)
}
p1_1 <- plot_strains_trait(gc.prm, lag, "lag time (hr)") + theme(axis.title.x = element_blank())
p1_2 <- plot_strains_trait(gc.prm, r, expression(growth~rate(h^-1))) + theme(strip.text = element_blank(), axis.title.x = element_blank())
p1_3 <- plot_strains_trait(gc.prm, maxOD, expression(paste("max", "[", OD[600], "]"))) + theme(strip.text = element_blank())

# By site
plot_boxplot_pair <- function (tb, ytrait, ylab = "") {
    tb %>%
        ggplot() +
        geom_rect(data = tibble(strain_site_group = c("H", "L")), aes(fill = strain_site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_boxplot(aes(x = strain_site_group, y = {{ytrait}}), fill = "white", outlier.size = -1, color = "black") +
        geom_point(aes(x = strain_site_group, y = {{ytrait}}, group = strain, color = strain), shape = 21, size = 2, stroke = 1, fill = NA,
                   position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
        scale_color_manual(values = rep("black", 100)) +
        scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
        #scale_x_discrete(label = c("high elevation", "low elevation")) +
        facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high elevation", L="low elevation"))) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
            strip.background = element_rect(color = NA, fill = NA),
            strip.text = element_text(size = 10, color = "black"),
            axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0,5,0,0), "mm")
        ) +
        guides(color = "none") +
        labs(x = "", y = ylab)

}
p2_1 <- plot_boxplot_pair(gc.prm, lag, "lag time (hr)") + theme(axis.title.x = element_blank())
p2_2 <- plot_boxplot_pair(gc.prm, r, expression(growth~rate(h^-1))) + theme(strip.text = element_blank(), axis.title.x = element_blank())
p2_3 <- plot_boxplot_pair(gc.prm, maxOD, expression(paste("max", "[", OD[600], "]"))) + theme(strip.text = element_blank(), axis.title.x = element_blank())

p1 <- plot_grid(p1_1, p2_1, ncol = 2, axis = "tbrl", align = "hv", rel_widths = c(2,1), labels = c("A", "B"), scale = .95)
p2 <- plot_grid(p1_2, p2_2, ncol = 2, axis = "tbrl", align = "hv", rel_widths = c(2,1), scale = .95)
p3 <- plot_grid(p1_3, p2_3, ncol = 2, axis = "tbrl", align = "hv", rel_widths = c(2,1), scale = .95)
p <- plot_grid(p1, p2, p3, ncol = 1, axis = "tbrl", align = "hv") + paint_white_background()


ggsave(here::here("plots/FigS8.png"), p, width = 10, height = 8)


## Does the rhizobia strain have effect on any growth trait?
mod <- lmer(lag ~ strain + (1|strain_site_group) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group has effect on lag time
# Response: lag
#               Chisq Df Pr(>Chisq)
# (Intercept) 130.457  1  < 2.2e-16 ***
#     strain    90.501 18  1.174e-11 ***
mod <- lmer(r ~ strain + (1|strain_site_group) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group has effect on r
# Response: r
#               Chisq Df Pr(>Chisq)
# (Intercept)  28.932  1  7.498e-08 ***
# strain      120.438 18  < 2.2e-16 ***
mod <- lmer(maxOD ~ strain + (1|strain_site_group) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group has effect on maxOD
# Response: maxOD
#           Chisq Df Pr(>Chisq)
# (Intercept)  53.848  1  2.166e-13 ***
#     strain  414.682 18  < 2.2e-16 ***

## Does the rhizobia sites have effect on any growth trait?
mod <- lmer(lag ~ strain_site_group + (1|strain) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group does not have effect on lag time
# Response: lag
#           Chisq Df Pr(>Chisq)
# (Intercept)       292.6038  1     <2e-16 ***
# strain_site_group   0.4997  1     0.4796
mod <- lmer(r ~ strain_site_group + (1|strain) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group does not have effect on r
# Response: r
#           Chisq Df Pr(>Chisq)
# (Intercept)       45.5158  1  1.514e-11 ***
# strain_site_group  1.7946  1     0.1804
mod <- lmer(maxOD ~ strain_site_group + (1|strain) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group does not have effect on maxOD
# Response: maxOD
# Chisq Df Pr(>Chisq)
# (Intercept)       36.7177  1  1.365e-09 ***
#  strain_site_group  1.8372  1     0.1753







