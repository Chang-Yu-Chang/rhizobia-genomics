#

library(tidyverse)
library(flextable) # for making table
source(here::here("analysis/00-metadata.R"))

features <- tibble(
    `Feature type` = c(rep("host yield", 2), rep("nodule", 2), rep("root architecture", 11)),
    `Feature` = traits
) %>%
    mutate(`Description` = c(
        "above-ground dry weight (mg)",
        "number of nodule",
        "root dry weight (mg)",
        "nodule dry weight (mg)",
        "number of root tips",
        "number of branching points",
        "total root length (px)",
        "branching frequency per px (1/px)",
        "area size (px^2)",
        "average diameter (px)",
        "median diameter (px)",
        "maximum diameter (px)",
        "estimated perimeter (px)",
        "estimated volume (px^3)",
        "estimated surface area (px^2)"
    )) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())

ft <- features %>%
    flextable() %>%
    merge_at(i = 1:2, j = 2, part = "body") %>%
    merge_at(i = 3:4, j = 2, part = "body") %>%
    merge_at(i = 5:15, j = 2, part = "body") %>%
    width(j = 2, width = 2) %>%
    width(j = 3, width = 2) %>%
    width(j = 4, width = 4) %>%
    #valign(j = 2, valign = "top", part = "all") %>%
    hline(i = 2, j = NULL, border = NULL, part = "body") %>%
    hline(i = 4, j = NULL, border = NULL, part = "body") %>%
    hline(i = 15, j = NULL, border = NULL, part = "body") %>%
    fix_border_issues()

save_as_image(ft, here::here("plots/TableS1.png"), webshot = "webshot2")


