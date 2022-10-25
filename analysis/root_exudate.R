#' Analysis of the root exudate data from Wilson et al 2015 where they have
#' they conducted transcriptomics and metabolomics of root border cells (BC) versus root tip (RT).
#' Here I only need the metabolomic data of BCs

library(tidyverse)

gcms <- read_csv(here::here("data/raw/root_exudate/253054TABLE_SIII_Watson_et_al_Border_Cells_Jan15.csv"), show_col_types = F)

# 00. clean up column names ----
#' A few concerns
#' 1. I am not sure if ...3 is the batch ID or means what. Given that the same CS
#' appears multiple times in one trial, I guess it should be
gcms <- gcms %>%
    select(...3, ...4, starts_with("BC"), starts_with("RT"), `avg BC`, `avg RT`) %>%
    rename(Batch = ...3, CS = ...4) %>%
    select(Batch, CS, "BC1...15", "BC2...16", "BC3...17", "BC4...18", "RT1...20", "RT2...21", "RT3...22", "RT4...23", "BC/RT", `avg BC`, `avg RT`) %>%
    rename_all(function(x){str_replace(x, "\\.\\.\\.\\d+$", "")}) %>%
    slice(-1) %>%
    filter(!is.na(CS), CS != "Unknown", CS != 1) %>%
    filter(CS != "NA)") %>%
    filter(!is.na(`BC/RT`)) %>%
    mutate(`BC/RT` = as.numeric(`BC/RT`)) %>%
    #distinct(Batch, CS, .keep_all = T) %>%
    arrange(CS, Batch)


gcms_clean <- gcms %>%
    select(Batch, CS, `BC/RT`, `avg BC`, `avg RT`) %>%
    arrange(desc(`avg BC`)) %>%
    # Remove amino acids
    filter(!str_detect(CS, "ine")) %>%
    filter(!str_detect(CS, "glutami")) %>%
    filter(!str_detect(CS, "Aspartic")) %>%
    filter(!str_detect(CS, "Ethylbis")) %>%
    # Remove inorganic acids
    filter(!str_detect(CS, "Phosphoric Acid")) %>%
    filter(!str_detect(CS, "Boric Acid")) %>%
    # Others
    filter(!str_detect(CS, "Urea")) %>%
    filter(!str_detect(CS, "Butanoic acid")) %>%
    filter(!str_detect(CS, "Myo-Inositol"))

coeff = 10
p1 <- gcms_clean %>%
    arrange(desc(`avg BC`)) %>%
    distinct(Batch, CS, .keep_all = T) %>%
    slice(1:10) %>%
    mutate(CS = factor(CS, CS)) %>%
    ggplot() +
    geom_col(aes(x = CS, y = `avg BC`), color = 1, fill = NA) +
    geom_point(aes(x = CS, y = `BC/RT` * coeff), size = 3, shape = 3) +
    geom_hline(yintercept = 1*coeff, linetype = 2, color = "red") +
    scale_y_continuous(name = "Average BC", sec.axis = sec_axis(~./coeff, name = "BC/RT", breaks = 1:10)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p2 <- gcms_clean %>%
    arrange(desc(`BC/RT`)) %>%
    distinct(Batch, CS, .keep_all = T) %>%
    slice(1:10) %>%
    mutate(CS = factor(CS, CS)) %>%
    ggplot() +
    geom_col(aes(x = CS, y = `avg BC`), color = 1, fill = NA) +
    geom_point(aes(x = CS, y = `BC/RT` * coeff), size = 3, shape = 3) +
    geom_hline(yintercept = 1*coeff, linetype = 2, color = "red") +
    scale_y_continuous(name = "Average BC", sec.axis = sec_axis(~./coeff, name = "BC/RT", breaks = 1:10)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p1
ggsave(here::here("plots/root_exudate_GCMS.png"), p, width = 5, height = 3)









