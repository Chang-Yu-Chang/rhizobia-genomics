#' This script plots growth curve traits

library(tidyverse)
library(janitor)
library(cowplot)
library(vegan) # for permanova
source(here::here("analysis/00-metadata.R"))

gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L")) %>%
    mutate(strain_site = str_sub(strain, 1, 2), strain_site_group = str_sub(strain, 1, 1))


# Clean up data
treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))

subset_ensifer <- function(tb) {
    tb %>%
        left_join(select(isolates_RDP, strain, Genus)) %>%
        drop_na()
}

gc.prm <- gc.prm %>% subset_ensifer()

# Panel A: cartoon for methods ----
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + draw_text("placeholder for\ncartoon")


# Panel B: growth rate at 30C ----
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
            #panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
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

p2 <- plot_boxplot_pair(gc.prm, r, expression(growth~rate(h^-1))) + theme(axis.title.x = element_blank())


# Panel C: PCA for growth traits ----
tt <- gc.prm %>%
    select(well, strain, strain_site_group, all_of(c("r", "lag", "maxOD"))) %>%
    drop_na()

pcobj <- tt %>%
    select(-well, -strain, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)
df <- as_tibble(predict(pcobj)[,1:2])

p3 <- df %>%
    bind_cols(select(tt, strain, strain_site_group)) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = strain_site_group)) +
    geom_polygon(stat = "ellipse", aes(x = PC1, y = PC2, color = strain_site_group), fill = NA, alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    scale_shape_manual(values = c(H = 16, L = 17), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = c(0.2, 0.1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = paste0("PC1 (", round(summary(pcobj)$importance[2,1]* 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pcobj)$importance[2,2]* 100, 1), "%)"))


## Stats
## Explanation of PC1 + PC2
summ <- summary(pcobj)
sum(summ$importance[2,1:2]) # 0.86127

# PERMANOVA test
Y <- gc.prm %>% select(c("r", "lag", "maxOD"))
set.seed(1)
adonis2(Y ~ strain, data = gc.prm, permutations = 999)
# adonis2(formula = Y ~ strain, data = gc.prm, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)
# strain   18  0.47521 0.67075 6.4513  0.001 ***
# Residual 57  0.23326 0.32925
# Total    75  0.70848 1.00000
adonis2(Y ~ strain_site_group, data = gc.prm, permutations = 999)
# adonis2(formula = Y ~ strain_site_group, data = gc.prm, permutations = 999)
#                    Df SumOfSqs      R2      F Pr(>F)
# strain_site_group  1  0.01154 0.01629 1.2257  0.263
# Residual          74  0.69693 0.98371
# Total             75  0.70848 1.00000


p <- plot_grid(p1, p2, p3, nrow = 1, axis = "tblr", align = "h", labels = LETTERS[1:3], scale = 0.95, rel_widths = c(1,1.5,1.5)) + paint_white_background()

ggsave(here::here("plots/Fig2.png"), p, width = 10, height = 4)














