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
plants <- read_csv(paste0(folder_data, "temp/23-plants.csv"))
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
#' Plot raw data and sem

compute_trait_mean <- function (plants_long, tra = "dry_weight_mg", pop = "VA") {
    pl <- plants_long %>%
        filter(trait == tra) %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, exp_id, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = exp_id, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = exp_id, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = exp_id, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
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
}

t1 <- compute_trait_mean(plants_long, tra = "dry_weight_mg", pop = "VA")
t2 <- compute_trait_mean(plants_long, tra = "dry_weight_mg", pop = "PA")
p1 <- plot_dots(t1$pl, t1$plm)
p2 <- plot_dots(t2$pl, t2$plm)
p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-04-trait_strain.png"), p, width = 5, height = 3)

# Stats
plants_test <- filter(plants, population == "VA")
## Does rhizobia strain have effect on shoot biomass in VA?
mod <- lmer(dry_weight_mg ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yup
## Does rhizobia strain have effect on nodule number in VA?
mod <- lmer(nodule_number ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yup
## Does rhizobia strain have effect on root biomass in VA?
mod <- lmer(root_weight_mg ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # nope

plants_test <- filter(plants, population == "PA")
## Does rhizobia strain have effect on shoot biomass in PA?
mod <- lmer(dry_weight_mg ~ exp_id + (1|site_group) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes
## Does rhizobia strain have effect on nodule number in VA?
mod <- lmer(nodule_number ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes
## Does rhizobia strain have effect on root biomass in VA?
mod <- lmer(root_weight_mg ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # nope


# 5. Plot the symbiosis traits comparing the two populations
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
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.2) +
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
}

t1 <- compute_trait_mean2(plants_long,tra = "dry_weight_mg", pop = "VA")
t2 <- compute_trait_mean2(plants_long,tra = "dry_weight_mg", pop = "PA")
p1 <- plot_dots2(t1$pl, t1$plm)
p2 <- plot_dots2(t2$pl, t2$plm)

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-05-trait_population.png"), p, width = 5, height = 3)

# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "VA")
mod <- lmer(dry_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # no
mod <- lmer(nodule_number ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # no
mod <- lmer(root_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # no

# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "PA")
mod <- lmer(dry_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes
mod <- lmer(nodule_number ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes
mod <- lmer(root_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes

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


# 7. Plot all three traits in dot plots by strain
compute_trait_mean3 <- function (plants_long, pop = "VA") {
    pl <- plants_long %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, exp_id, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots3 <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = exp_id, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = exp_id, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = exp_id, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
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
}

t1 <- compute_trait_mean3(plants_long, pop = "VA")
t2 <- compute_trait_mean3(plants_long, pop = "PA")
p1 <- plot_dots3(t1$pl, t1$plm)
p2 <- plot_dots3(t2$pl, t2$plm)

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-07-trait_strain.png"), p, width = 8, height = 9)


# 6. Plot all three traits in violin plot
compute_trait_mean4 <- function (plants_long, pop = "VA") {
    pl <- plants_long %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots4 <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.2) +
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
}

t1 <- compute_trait_mean4(plants_long, "VA")
t2 <- compute_trait_mean4(plants_long, "PA")
p1 <- plot_dots4(t1$pl, t1$plm)
p2 <- plot_dots4(t2$pl, t2$plm)



p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/23a-08-trait_population.png"), p, width = 8, height = 9)

