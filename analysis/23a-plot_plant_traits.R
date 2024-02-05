#' This script plots the trait data of plant inoculation experiments

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
library(ellipse) # for calculating the ellipse
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(broom.mixed) # for tidy up lme4
library(factoextra) # for plotting pca eclipse
source(here::here("analysis/00-metadata.R"))

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
options(contrasts=c("contr.sum", "contr.poly"))

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
sysplants <- read_csv(paste0(folder_data, "temp/23-plants.csv"))
plants_long <- read_csv(paste0(folder_data, "temp/23-plants_long.csv"))
plants_wide <- read_csv(paste0(folder_data, "temp/23-plants_wide.csv"))

site_group_colors <- c(`high elevation` = "#0C6291", `low elevation` = "#BF4342", 
                 `suburban` = "#0C6291", `urban` = "#BF4342", control = "grey")

# 1. plot all plant traits
length(unique(plants_long$id)) # 231 plants
length(unique(plants_long$exp_id)) # 14 rhizobia + 1 control
plants_long %>%
    arrange(site_group) %>%
    distinct(id, site_group, exp_id) %>%
    group_by(site_group, exp_id) %>%
    count()

p <- plants_long %>%
    ggplot() +
    geom_boxplot(aes(x = exp_id, y = value, color = site_group), outlier.shape = NA) +
    geom_jitter(aes(x = exp_id, y = value, color = site_group), width = 0.2, height = 0, shape = 21) +
    facet_grid(trait~site_group, scales = "free", space = "free_x") +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    guides(color = "none") +
    labs()

ggsave(paste0(folder_data, "temp/23a-01-plant_traits.png"), p, width = 10, height = 10)

# 2. Plot the nodule vs plant biomass
p <- plants_wide %>%
    filter(site_group != "control") %>%
    ggplot() +
    geom_smooth(aes(x = nodule_number, y = dry_weight_mg), method = "lm")+
    geom_point(aes(x = nodule_number, y = dry_weight_mg), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
    #scale_color_manual(values = brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)]) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    coord_equal() +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey95")
    ) +
    guides() +
    labs(x = "# of nodules", y = "shoot biomass (mg)")

ggsave(paste0(folder_data, "temp/23a-02-nodule_vs_shoot.png"), p, width = 5, height = 4)

# Stats
lm(dry_weight_mg ~ nodule_number, data = plants_wide) %>%
    broom::tidy()

# 3. Plot the pairwise comparison of the three symbiosis traits
plot_dots <- function(tb_wide, trait1, trait2) {
    tb_wide %>%
        ggplot() +
        geom_smooth(aes(x = {{trait1}}, y = {{trait2}}), method = "lm")+
        geom_point(aes(x = {{trait1}}, y = {{trait2}}), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
        #scale_color_manual(values = brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)]) +
        #scale_x_continuous(breaks = seq(0, 50, 10)) +
        #facet_wrap(~site_group, nrow = 2) +
        #coord_equal() +
        theme_classic() +
        theme(
            panel.grid.major = element_line(color = "grey95")
        ) +
        guides() 
    #    labs(x = "# of nodules", y = "shoot biomass (mg)")
}
p1 <- plants_wide %>%
    filter(site_group != "control") %>%
    plot_dots(nodule_number, dry_weight_mg) +
    labs(x = "# of nodules", y = "shoot biomass (mg)")
p2 <- plants_wide %>%
    filter(site_group != "control") %>%
    plot_dots(nodule_number, root_weight_mg) +
    labs(x = "# of nodules", y = "root biomass (mg)")
p3 <- plants_wide %>%
    filter(site_group != "control") %>%
    plot_dots(root_weight_mg, dry_weight_mg) +
    labs(x = "root biomass (mg)", y = "shoot biomass (mg)")
p <- plot_grid(p1, p2, p3,
    nrow = 1, align = "h", axis = "tb", scale = 0.95, labels = LETTERS[1:3]
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(paste0(folder_data, "temp/23a-03-symbiosis_trait_pairs.png"), p, width = 9, height = 3)

# Are the symbiosis traits mutually correlated?
cor.test(plants_wide$nodule_number, plants_wide$dry_weight_mg) %>% broom::tidy()
cor.test(plants_wide$nodule_number, plants_wide$root_weight_mg) %>% broom::tidy()
cor.test(plants_wide$root_weight_mg, plants_wide$dry_weight_mg) %>% broom::tidy()

# 4. Plot the symbiosis traits comparing the strains
p1 <- plants_long %>% 
    filter(trait == "dry_weight_mg") %>%
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = exp_id, y = value), fill = "white") +
    geom_boxplot(aes(x = exp_id, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "shoot biomass (mg)")
p2 <- plants_long %>% 
    filter(trait == "dry_weight_mg") %>%
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = exp_id, y = value), fill = "white") +
    geom_boxplot(aes(x = exp_id, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "shoot biomass (mg)")

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-04-trait_strain.png"), p, width = 5, height = 3)

# Does rhizobia strain have effect on shoot biomass in VA?
plants_test <- filter(plants, population == "VA")
mod <- lmer(dry_weight_mg ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # Difference In PA

# Does rhizobia strain have effect on shoot biomass in PA?
plants_test <- filter(plants, population == "PA")
mod <- lmer(dry_weight_mg ~ exp_id + (1|site_group) + (1|waterblock), data = plants_test)
Anova(mod, type = 3)  # Differnece berween urban and suburban


# 5. Plot the symbiosis traits comparing the two populations
p1 <- plants_long %>% 
    filter(trait == "dry_weight_mg") %>%
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = site_group, y = value), fill = "white") +
    geom_boxplot(aes(x = site_group, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank()    
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "shoot biomass (mg)")

p2 <- plants_long %>% 
    filter(trait == "dry_weight_mg") %>%
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = site_group, y = value), fill = "white") +
    geom_boxplot(aes(x = site_group, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank()    
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "shoot biomass (mg)")

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-05-trait_population.png"), p, width = 5, height = 3)

# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "VA")
mod <- lmer(dry_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # No difference between high and low elevation
mod <- lmer(nodule_number ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # No difference between high and low elevation
mod <- lmer(root_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # No difference between high and low elevation


# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "PA")
mod <- lmer(dry_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3)  # Differnece berween urban and suburban
mod <- lmer(nodule_number ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3)  # Differnece berween urban and suburban
mod <- lmer(root_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3)  # Differnece berween urban and suburban

# 6. Plot the PCA for the four sites
pw1 <- plants_wide %>%
    filter(population == "VA", exp_id != "control")
pcobj1 <- pw1 %>%
    select(dry_weight_mg, nodule_number, root_weight_mg) %>%
    prcomp(center = TRUE, scale. = TRUE)

p1 <- fviz_pca_ind(
    pcobj1,
    label = "none",
    habillage = pw1$site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "top",
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()

pw2 <- plants_wide %>%
    filter(population == "PA", exp_id != "control")
pcobj2 <- pw2 %>%
    select(dry_weight_mg, nodule_number, root_weight_mg) %>%
    prcomp(center = TRUE, scale. = TRUE)

p2 <- fviz_pca_ind(
    pcobj2,
    label = "none",
    habillage = pw2$site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "top",
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()


p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"), scale = 0.9) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/23a-06-trait_pca.png"), p, width = 9, height = 5)

# Does rhizobia sites have effect on PC1 of symbiosis in VA?
plants_test <- mutate(pw1, pc1 = pcobj1$x[,1]) %>% left_join(plants)
mod <- lmer(pc1 ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # No difference between high and low elevation

# Does rhizobia sites have effect on PC1 of symbiosis in VA?
plants_test <- mutate(pw2, pc1 = pcobj2$x[,1]) %>% left_join(plants)
mod <- lmer(pc1 ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3)  # No differnece between urban and suburban
#pcobj2$sdev^2 / sum(pcobj2$sdev^2)



# 7. Plot all three traits in violin plot
p1 <- plants_long %>% 
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = exp_id, y = value), fill = "white") +
    geom_boxplot(aes(x = exp_id, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_grid(trait~site_group, scales = "free") +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "")

p2 <- plants_long %>% 
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = exp_id, y = value), fill = "white") +
    geom_boxplot(aes(x = exp_id, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_grid(trait~site_group, scales = "free") +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "")

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-07-trait_strain.png"), p, width = 8, height = 9)



# 7. Plot all three traits in violin plot
p1 <- plants_long %>% 
    filter(exp_id != "control") %>%
    filter(population == "VA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = site_group, y = value), fill = "white") +
    geom_boxplot(aes(x = site_group, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_grid(trait~site_group, scales = "free") +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "")

p2 <- plants_long %>% 
    filter(exp_id != "control") %>%
    filter(population == "PA") %>%
    ggplot() +
    geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
    geom_violin(aes(x = site_group, y = value), fill = "white") +
    geom_boxplot(aes(x = site_group, y = value), width=0.1) +
    scale_fill_manual(values = site_group_colors) +
    facet_grid(trait~site_group, scales = "free") +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(x = " ", y = "")

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-09-trait_population.png"), p, width = 8, height = 9)


if (FALSE) {









# Stat for slope in regression
lm_x_y <- function (dat, x_var, y_var) {
    formula <- as.formula(paste(y_var, "~", x_var))
    lm(formula, data = dat) %>%
        broom::tidy() %>%
        filter(term == x_var) %>%
        select(term, estimate, p.value)
}

plants_wide %>%
    filter(site_group != "control") %>%
    nest(data = c(id, exp_id, dry_weight_mg, nodule_number, root_weight_mg)) %>%
    mutate(lm = map(data, ~lm_x_y(.x, "nodule_number", "dry_weight_mg"))) %>%
    select(site_group, lm) %>%
    unnest(lm)

# 4. plot the root traits 
p <- plants %>%
    #select(id, exp_id, site_group, all_of(traits2)) %>%
    #select(id, site_group, exp_id, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>%
    pivot_longer(cols = c(-id, -exp_id, -site_group), names_to = "trait") %>%
    ggplot() +
    geom_boxplot(aes(x = exp_id, y = value, color = site_group), outlier.shape = NA) +
    geom_jitter(aes(x = exp_id, y = value, color = site_group), width = 0.2, height = 0, shape = 21) +
    facet_grid(trait~site_group, scales = "free", space = "free_x") +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme_facets +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    guides(color = "none") +
    labs()

ggsave(paste0(folder_data, "temp/23a-04-root_traits.png"), p, width = 8, height = 40)

# 5. plot the trait correlation 
# Names
plants_traits <- plants %>%
    filter(site_group %in% c("high-elevation", "low-elevation")) %>%
    select(all_of(traits))
tb_traits <- tibble(trait = ordered(names(plants_traits), names(plants_traits))) %>%
    expand(trait, trait) %>%
    rename(trait1 = trait...1, trait2 = trait...2) %>%
    filter(trait1 < trait2) %>%
    arrange(trait1, trait2)

# Calculate cor
cor_matrix <- cor(plants_traits, use = "complete.obs")
cor_x_y <- function (dat) cor.test(unlist(dat[,1]), unlist(dat[,2]), method = "spearman", exact = F) %>% broom::tidy() %>% clean_names()
tb_cortest <- tb_traits %>%
    rowwise() %>%
    mutate(traits = list(select(plants_traits, trait1, trait2))) %>%
    ungroup() %>%
    # Compute cor test
    mutate(cortest = map(traits, ~cor_x_y(.x))) %>%
    select(-traits) %>%
    unnest(cortest) %>%
    select(trait1, trait2, cor=estimate, p_value, method, alternative) %>%
    mutate(asterisk = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        T ~ "n.s."
    ))

# tb_cor <- as_tibble(cor_matrix) %>%
#     mutate(trait1 = colnames(.)) %>%
#     pivot_longer(-trait1, names_to = "trait2", values_to = "cor") %>%
#     mutate(trait1 = ordered(trait1, names(plants_traits))) %>%
#     mutate(trait2 = ordered(trait2, names(plants_traits))) %>%
#     filter(trait1 < trait2) %>%
#     arrange(trait1, trait2)




# Calculate ellipse
tb_ell <- tb_traits %>%
    rowwise() %>%
    # Calculate ellipse
    mutate(cor_m = list(ellipse(cor_matrix, which = which(names(plants_traits) %in% c(trait1, trait2))))) %>%
    mutate(cor_tb = list(tibble(x = cor_m[,which(colnames(cor_m) == trait1)], y = cor_m[,which(colnames(cor_m) == trait2)]))) %>%
    select(-cor_m) %>%
    unnest(cor_tb)

tb_cortest_swap <- tibble(trait1 = tb_cortest$trait2, trait2 = tb_cortest$trait1, cor = tb_cortest$cor, asterisk = tb_cortest$asterisk)

#
p <- tb_ell %>%
    left_join(tb_cortest) %>%
    ggplot() +
    geom_polygon(aes(x = x, y = y, fill = cor), color = "black") +
    geom_text(data = tb_cortest_swap, x = 0, y = 0, aes(label = round(cor, 2))) +
    geom_text(data = tb_cortest_swap, x = 0, y = 0, aes(label = asterisk), vjust = -0.5) +
    facet_grid(trait2 ~ trait1, switch = "y", drop = F) +
    scale_fill_gradient2(low = "maroon", high = "steelblue", breaks = c(1,0,-1), limits = c(-1,1)) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        #panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_rect(color = "white", fill = "white")
    ) +
    guides() +
    labs()


ggsave(paste0(folder_data, "temp/23a-05-trait_cor.png"), p, width = 12, height = 10)













}
