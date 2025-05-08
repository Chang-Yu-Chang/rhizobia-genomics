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
gtw <- read_csv(paste0(folder_phenotypes, 'growth/gtw.csv')) %>% # Growth traits per well
    left_join(select(iso, exp_id, contig_species))

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
        p_em = map2(em, trait, plot_emmeans)
    )

p <- plot_grid(plotlist = tb_mod$p_em, ncol = 1)

ggsave(paste0(folder_data, "phenotypes/growth/emmeans.png"), p, width = 6, height = 10)

tb_mod$mod_tidied[[1]]
tb_mod$mod_tidied[[2]]
tb_mod$mod_tidied[[3]]

#check_model(tb_mod$mod[[3]])
