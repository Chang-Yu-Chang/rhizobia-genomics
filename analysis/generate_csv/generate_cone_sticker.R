#' This script reads the generate the tabular csv for printing stickers for cone containers

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# 1. Count the number of plants ----
#' The number of plants were counted the day before rhizobia inoculation
rhizobia_mid <- c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2")
rhizobia_HL <- c("H3M1R1", "L2M2R1")

plants_count <- tibble(
    Plant = c("S5M1", "S5M3", "S5M5", "S5M4", "S5M6", "S5M7",
              "H2M1", "H2M2", "H2M3", "H4M1", "H4M3",
              "L1M1", "L1M3", "L3M1", "L3M2"),
    Count = c(15, 26, 12, 32, 21, 10,
              10, 2, 3, 6, 5,
              6, 8, 2, 9)) %>%
    mutate(PlantSite = str_sub(Plant, 1, 1)) %>%
    arrange(Plant)


# Mid sites
plants_mid <- plants_count %>%
    filter(PlantSite == "S") %>%
    # Replicate the rows
    rowwise() %>%
    split.data.frame(f = .$Plant) %>%
    lapply(function(x) {x[rep(1, x$Count),]}) %>%
    bind_rows() %>%
    select(-Count) %>%
    # Assign rhizobia equally
    group_by(Plant) %>%
    mutate(Rhizobia = c(rep(rhizobia_mid, n() %/% length(rhizobia_mid)), rep(NA, n() %% length(rhizobia_mid)))) %>%
    arrange(Rhizobia, Plant) %>%
    mutate(RhizobiaSite = str_sub(Rhizobia, 1, 1)) %>%
    select(RhizobiaSite, Rhizobia, PlantSite, Plant) %>%
    ungroup()

# Mid sites
plants_HL <- plants_count %>%
    filter(PlantSite %in% c("H", "L")) %>%
    # Replicate the rows
    rowwise() %>%
    split.data.frame(f = .$Plant) %>%
    lapply(function(x) {x[rep(1, x$Count),]}) %>%
    bind_rows() %>%
    select(-Count) %>%
    # Assign rhizobia equally
    group_by(Plant) %>%
    mutate(Rhizobia = c(rep(rhizobia_HL, n() %/% length(rhizobia_HL)), rep(NA, n() %% length(rhizobia_HL)))) %>%
    arrange(Rhizobia, Plant) %>%
    mutate(RhizobiaSite = str_sub(Rhizobia, 1, 1)) %>%
    select(RhizobiaSite, Rhizobia, PlantSite, Plant) %>%
    ungroup()

# 2. Generate treatment ID and assign waterblock ----
#' In this expriment the plants are grouped into a unit of 6 as the blue water block can hold 6 plants
# Treatment ID
treatment_ID <- bind_rows(plants_mid, plants_HL) %>%
    group_by(Rhizobia, PlantSite) %>%
    distinct(Rhizobia, PlantSite) %>%
    arrange(Rhizobia, PlantSite) %>%
    ungroup() %>%
    filter(!is.na(Rhizobia)) %>%
    mutate(TreatmentID = 1:n()) %>%
    select(TreatmentID, everything())

treatments_temp <- bind_rows(plants_mid, plants_HL) %>%
    left_join(treatment_ID) %>%
    group_by(TreatmentID) %>%
    group_split

# Assign waterblockID
counter <- 1
temp <- rep(list(NA, length(treatments_temp)))
set.seed(1)
for (i in 1:length(treatments_temp)) {
    temp[[i]] <- treatments_temp[[i]] %>%
        slice(order(runif(n()))) %>%
        mutate(Waterblock = `[`(rep(counter:(counter + nrow(treatments_temp[[i]]) %/% 6), each = 6), 1:n()))
    counter <- counter + nrow(treatments_temp[[i]]) %/% 6 + 1
}

treatments <- temp %>%
    bind_rows() %>%
    mutate(ID = 1:n()) %>%
    ungroup() %>%
    mutate(ID = sprintf("%03d", 1:n())) %>%
    mutate(Label = paste0("ID", ID, "-", Rhizobia, "-WB", Waterblock, "-", Plant)) %>%
    select(ID, everything())

# Uncomment this line to generate the csv table
#write_csv(treatments, paste0(folder_data, "raw/rhizobia/03-label/treatments.csv"))

#' To use the printing labels in MS Words, save the csv file into a xlsx file
#' For how to print the labels, refer to the Google Docs from Corlett
#' "Create and print pot labels with Excel and Mail Merge in Word"

# 3. Manually assign the rest of the plants that are not included in the treatments.csv ----
#' There were a few plants manually assigned during the day of inoculation
#' Specifically, there are 17 of them
#' To include these plants in the data sheet, append the following table to treatments
treatments_unassigned <- tibble(
    ID = sprintf("%02d", 151:167),
    RhizobiaSite = c(rep("H",4), rep("L",3), "H", rep(NA, 7), "L", NA),
    Rhizobia = c("H3M1R1", rhizobia_mid, "H3M1R1", rep(NA, 7), "L2M2R1", NA),
    PlantSite = c("L", rep("S", 6), "H", rep("S", 7), "H", "S"),
    Plant = treatments$Plant[151:167],
    Waterblock = c(6,3,12,15,24,27,30,6, rep(31,4), rep(32,3), 9, 32),
    TreatmentID = c(3, 1,4,5, 8,9,10, 2, rep(0, 7), 6, 0)
) %>%
    mutate(Label = paste0("ID", ID, "-", Rhizobia, "-WB", Waterblock, "-", Plant))

treatments_assigned <- treatments %>%
    rows_update(treatments_unassigned, by = "ID")

# Uncomment this line to generate the csv table
#write_csv(treatments_assigned, paste0(folder_data, "raw/rhizobia/03-label/treatments_assigned.csv"))



# 4. Count sample size per treatment ----
treatments_assigned %>%
    group_by(TreatmentID) %>%
    count()


