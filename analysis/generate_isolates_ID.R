#' This scripts generate the table of isolates with glycerol stocks

library(tidyverse)

# 1. Generate local-urban isolate ID  ----

isolates <- tibble(
    IsolateID = c(
        # 20221017
        "L2M4R1", "L1M3R2", "L1M3R1", "L1M2R2", "L1M1R4", "L2M4R4", "L2M4R2", "L2M3R2", "L2M2R1", "L2M1R1",
        "H4M5R2", "H4M1R1", "H3M4R1", "H3M3R2", "H3M3R1", "H3M2R1", "H3M1R1", "H2M3R1", "H2M2R1", "H4M5R1", "H2M4R1", "H2M3R2", "H1M2R3", "H1M2R1", "H1M1R1",
        # 20221019
        "L1M1R1", "L1M1R2", "L2M5R2", "L3M4R1",
        "H1M3R1", "H1M4R1", "H1M2R2", "H2M2R2", "H3M3R4", "H3M4R2",
        "L4M2R2", "L4M3R3", "L4M4R1", "L4M4R3", "L4M5R1",
        "L3M2R1", "L3M2R2", "L3M3R1", "L3M5R1", "L3M5R2",
        "L3M1R1", "L4M5R4", "L4M7R1", "L4M7R2", "L4M8R1",
        "L4M5R2", "L3M6R2", "L3M6R1", "L3M4R3", "L3M4R2", "L2M2R2", "L1M4R2", "L1M2R4", "L1M2R3",
        "H4M1R3",
        # 20221020
        "H4M2R2", "H4M4R2",
        "L1M2R1", "L1M4R1", "L4M1R1", "L4M1R3", "L4M5R3"
    )
) %>%
    mutate(Site = str_sub(IsolateID, 1, 1), SiteID = str_sub(IsolateID, 2, 2), PlantID = str_sub(IsolateID, 1, 4)) %>%
    select(Site, SiteID, PlantID, IsolateID)

isolates %>%
    distinct(IsolateID)


isolates %>%
    arrange(IsolateID) %>%
    group_by(PlantID) %>%
    count(name = "Count")

# Number of isolates at H and L
isolates %>%
    group_by(Site) %>%
    count(name = "Count")

# Number per sites
isolates %>%
    group_by(Site, SiteID) %>%
    count(name = "Count")


isolates <- isolates %>%
    arrange(IsolateID) %>%
    mutate(` `= 1:n()) %>%
    select(` `, everything())


write_csv(isolates, file = here::here("data/raw/isolates_ID.csv"))


#
isolate_cryo_L <- isolates %>%
    filter(Site == "L") %>%
    arrange(IsolateID) %>%
    pull(IsolateID) %>%
    c(rep("", 9*((length(.)%/%9)+1)-length(.))) %>%
    matrix(ncol = 9, byrow = T) %>%
    as_tibble() %>%
    set_names(1:9)

isolate_cryo_H <- isolates %>%
    filter(Site == "H") %>%
    arrange(IsolateID) %>%
    pull(IsolateID) %>%
    c(rep("", 9*((length(.)%/%9)+1)-length(.))) %>%
    matrix(ncol = 9, byrow = T) %>%
    as_tibble() %>%
    set_names(1:9) %>%
    bind_rows(isolate_cryo_L[5,])

padding <- matrix("", ncol = 9) %>% as_tibble %>% set_names(1:9)
isolate_cryo_L <- bind_cols(
    tibble(` ` = 1:8),
    bind_rows(
        isolate_cryo_L[1,],
        padding,
        isolate_cryo_L[2,],
        padding,
        isolate_cryo_L[3,],
        padding,
        isolate_cryo_L[4,],
        padding
    )
)
isolate_cryo_H <- bind_cols(
    tibble(` ` = 1:8),
    bind_rows(
        isolate_cryo_H[1,],
        padding,
        isolate_cryo_H[2,],
        padding,
        isolate_cryo_H[3,],
        padding,
        isolate_cryo_H[4,],
        padding
    )
)
write_csv(isolate_cryo_L, here::here("data/raw/isolates_cryo_L.csv"))
write_csv(isolate_cryo_H, here::here("data/raw/isolates_cryo_H.csv"))



# 2. Combine ID for my isolates and Terrence's isolates ----
isolates_urban <- read_csv(here::here("data/raw/rhizobia/isolates_urban.csv"), show_col_types = F) %>%
    filter(!is.na(StrainID))

isolates_for_seq <- tibble(ExpID = c(isolates$IsolateID, isolates_urban$StrainID)) %>%
    mutate(`Sample Name` = sprintf("%03.0f", 1:n()),
           `# of picks per container` = 1,
           Note = paste0("bag ", c(rep(1:5, each = 16), rep(6, 11))))

write_csv(isolates_for_seq, here::here("data/raw/rhizobia/isolates_for_seq.csv"))


c("Sample Name", "# of picks per container", "Note")




