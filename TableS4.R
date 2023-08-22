#' This script makes the table of the medium recipe used in the current study

library(tidyverse)
library(flextable)
library(ftExtra) # for formatting the subscripts
library(gdtools) # for font


# Stock          | Amount     | Concentration
# ---------------|------------|---------------
# Tryptone       | 5 g        |  5 g/L
# Yeast extract  | 3 g        |  3 g/L
# CaCl~2~        | 1.47 g     |  1.47 g/L
# Agar           | 18 g       |  18 g/L
# ddH~2~O        | 1 L        |
#
tt <- tibble(
    Stock = c(
        "Tryptone",
        "Yeast extract",
        "CaCl~2~",
        "Agar"
    ),
    `Amount (g)` = c(5, 3, 1.47, 18),
    `Final conc. (g/L)` = c(5, 3, 1.47, 18)
)

ft <- tt %>%
    mutate(across(everything(), as.character)) %>%
    flextable() %>%
    colformat_md() %>%
    width(j = 1:2, width = 1.5) %>%
    width(j = 3, width = 2) %>%
    footnote(i = 4, j = 1, ref_symbols = "a",
             value = as_paragraph("For liquid medium, neglect the agar"))

save_as_image(ft, here::here("plots/TableS4.png"), webshot = "webshot2")







