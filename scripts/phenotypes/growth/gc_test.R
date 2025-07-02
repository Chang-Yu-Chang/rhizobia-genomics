#'

library(tidyverse)
library(cowplot)
library(lme4)
library(car) # for anova
library(emmeans) # for post hoc
library(performance) # for evaluating model performance
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    mutate(contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gtw <- read_csv(paste0(folder_phenotypes, 'growth/gtw.csv')) %>% # Growth traits per well
    left_join(select(iso, exp_id, contig_species))
dml <- read_csv(paste0(folder_phenotypes, "sites/dml.csv")) # daily max t at sampling sites
sites  <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>%
    # Use sites where the isolates were from
    filter(site %in% isolates$site)
plot_emmeans <- function (em, trait) {
    plot(em) +
        facet_grid(~temperature, scales = "free_x") +
        coord_flip() +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(title = trait)
}


# By speceis ----
tb_mod <- tibble(
    trait = c("r", "lag", "maxOD"),
    formula = c(
        "r ~ contig_species*temperature + (1|well)",
        "lag ~ contig_species*temperature + (1|well)",
        "maxOD ~ contig_species*temperature + (1|well)"
    )
) %>%
    mutate(
        mod = map(formula, ~lmer(as.formula(.x), data = gtw)),
        mod_tidied = map(mod, ~Anova(.x, type = 3)),
        em = map(mod, ~emmeans(.x, ~contig_species+temperature)),
        p_em = map2(em, trait, plot_emmeans),
        empair = map(em, ~broom::tidy(pairs(.x)))
    )

p <- plot_grid(plotlist = tb_mod$p_em, ncol = 1)

ggsave(paste0(folder_data, "phenotypes/growth/emmeans.png"), p, width = 6, height = 10)

#
tb_mod$mod_tidied[[1]]
tb_mod$mod_tidied[[2]]
tb_mod$mod_tidied[[3]]

tb_mod$empair[[1]] %>%
    #filter(contrast == "S. meliloti 25c - S. medicae 25c")
    #filter(contrast == "S. meliloti 30c - S. medicae 30c")
    #filter(contrast == "S. meliloti 35c - S. medicae 35c")
    filter(contrast == "S. meliloti 40c - S. medicae 40c")

tb_mod$empair[[3]] %>%
    filter(contrast == "S. meliloti 35c - S. canadensis 35c")


tb_mod$empair[[3]] %>%
    #filter(contrast == "S. meliloti 25c - S. medicae 25c")
    #filter(contrast == "S. meliloti 30c - S. medicae 30c")
    #filter(contrast == "S. meliloti 35c - S. medicae 35c")
    filter(contrast == "S. meliloti 40c - S. medicae 40c")

# by habitat tmax and site ----
dmli <- dml %>%
    group_by(population, site) %>%
    summarize(tmax_max = max(tmax_deg_c), tmin_max = max(tmin_deg_c))

tb <- gtw %>%
    select(exp_id, well, r, lag, maxOD, temperature) %>%
    left_join(select(iso, exp_id, contig_species)) %>%
    left_join(select(isolates, exp_id, population, site)) %>%
    left_join(dmli)
    #filter(contig_species %in% c("S. medicae", "S. meliloti"))

# tmax
mod <- lmer(r ~ temperature*tmax_max + contig_species + (1|exp_id:well) + (1|site), data = tb)
Anova(mod, type = 3)

mod <- lmer(lag ~ temperature*tmax_max + contig_species + (1|exp_id:well) + (1|site), data = tb)
Anova(mod, type = 3)

mod <- lmer(maxOD ~ temperature*tmax_max + contig_species  + (1|exp_id:well) + (1|site), data = tb)
Anova(mod, type = 3)

# tmin
mod <- lmer(r ~ temperature*tmin_max + contig_species + (1|exp_id:well) + (1|site), data = tb)
Anova(mod, type = 3)

mod <- lmer(lag ~ temperature*tmin_max + contig_species + (1|exp_id:well) + (1|site), data = tb)
Anova(mod, type = 3)
emmeans(mod, ~temperature) %>% pairs

mod <- lmer(maxOD ~ temperature*tmin_max + contig_species  + (1|exp_id:well) + (1|site), data = tb)
Anova(mod, type = 3)

# by habitat tmax, within region ----
tb_mod <- tibble(
    trait = rep(c("r", "lag", "maxOD"), each = 2),
    pop = rep(c("VA", "PA"), 3),
    formula = rep(c(
        "r ~ temperature*tmax_max + contig_species  + (1|exp_id:well) + (1|site)",
        "lag ~ temperature*tmax_max + contig_species  + (1|exp_id:well) + (1|site)",
        "maxOD ~ temperature*tmax_max + contig_species  + (1|exp_id:well) + (1|site)"
    ), each = 2)
) %>%
    mutate(
        mod = map2(formula, pop, ~lmer(as.formula(.x), data = filter(tb, population == .y))),
        mod_tidied = map(mod, ~Anova(.x, type = 3))
    )

# tmax
tb_mod$mod_tidied[[1]]
tb_mod$mod_tidied[[2]]
tb_mod$mod_tidied[[3]]
tb_mod$mod_tidied[[4]]
tb_mod$mod_tidied[[5]]
tb_mod$mod_tidied[[6]]
