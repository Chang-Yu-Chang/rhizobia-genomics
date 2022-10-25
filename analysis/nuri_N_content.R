#' Nuri's medicago data under fertilizer with different N content

library(tidyverse)

df <- read_csv(here::here("data/raw/nuri_N_content/nurirunfullfinal.csv"), show_col_types = F)
onetr <- read_csv(here::here("data/raw/nuri_N_content/1randomtechrepNRfinal.csv"), show_col_types = F)

# 0. clean up the data ----
onetr <- onetr %>%
    mutate(Total_biomass = Root_biomass + Shoot_biomass,
           scaledN = scale(N_percent),
           scaleNsq = scaledN^2,
           scaled_biomass = scale(Total_biomass)) %>%
    filter(!is.na(Number)) %>%
    filter(!is.na(Shoot_biomass)) %>%
    filter(!is.na(Root_biomass)) %>%
    filter(Mortality_End == "1")

# 0.1 Exclude the R- plants that formed nodules ----
#Making a dataframe for nodule analysis: excluding R- plants that formed nodules
#and also excluding #209, #210 (A145 originally) because they formed a lot of pink nodules
onetr <- onetr %>%
    filter(Rhizobia != "R-", Rstate == "R+") %>%
    filter(Number != 209 & Number != 210) %>%
    select(Rhizobia, Number, Total_nod)

# 0.2 Include R- plants that did form a few nodules ----
#this dataframe should be when nodule counts/types matter
Rplus1 <- onetr %>%
    filter(Rhizobia == "R-" & Total_nod >= 1) %>%
    mutate(Rhizobia = "Rplus")
#Plants 209, 210 (A145 originally) should be moved into this undefined Rplus category as well
#Because they formed a high number of pink nodules: assume contamination with a different R strain.
Rplus2 <- onetr %>%
    filter(Number %in% 209:210) %>%
    mutate(Rhizobia = "Rplus")

Rplusmerged <- bind_rows(onetr, Rplus1, Rplus2)

# 1. nodule analysis ----

