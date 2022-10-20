library(tidyverse)

seeds <- tibble(
    PlantID = c(
        "H1M1", "H1M3", "H1M4",
        "H2M1", "H2M2", "H2M3", "H2M4",
        "H3M1", "H3M2", "H3M3", "H3M4",
        "H4M1", "H4M2", "H4M3", "H4M4", "H4M5",
        "L1M1", "L1M3",
        "L2M6",
        "L3M1", "L3M2",
        "L3M4", "L3M7", "L4M8"
    ),

    SeedCount = c(
        36,23,4,
        55,86,33,44,
        55,72,14,3,
        49,22,86,18,17,
        10,23,
        6,
        11,20,
        7,13,8
    )
) %>%
    mutate(Site = str_sub(PlantID, 1, 1), SiteID = str_sub(PlantID, 2, 2)) %>%
    select(Site, SiteID, PlantID, SeedCount)

seeds %>%
    group_by(Site, SiteID) %>%
    summarize(SeedCount = sum(SeedCount))

seeds %>%
    group_by(Site, SiteID) %>%
    count()

