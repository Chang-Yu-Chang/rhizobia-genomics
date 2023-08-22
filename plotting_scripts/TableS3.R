#' This script makes the table of the medium recipe used in the current study

library(tidyverse)
library(flextable)
library(ftExtra) # for formatting the subscripts
library(gdtools) # for font

tt <- tibble(
    ` ` = c(rep("macronutrient", 5), rep("micronutrient", 5)),
    Stock = c(
        "MgSO~4~ 7H~2~O",
        "KH~2~PO~4~",
        "Na~2~HPO~4~ 2H~2~O",
        "CaCl~2~ 2H~2~O",
        "Fe-EDTA",

        "MnSO4",
        "CuSO4",
        "0.1M ZnSO~4~",
        "H~3~BO~3~",
        "Na~2~MoO~4~"
    ),
    `Molar mass (g/mole)` = c(
        120, 136, 142, 111, NA,
        151, 160, 161, 62, 206
    ),
    Ratio = c(
        120/(120+7*18), 1, 142/(142+2*18), 111/(111+2*18), NA,
        1,1,1,1,1
    ) %>% round(3),
    `Amount (g)` = c(
        30.8, 23.825, 17.8, 36.753, NA,
        0.01, 0.01, 619.4, 0.01, 0.01
    ),
    `Water (L)` = c(
        rep(0.25, 4), NA,
        rep(0.1, 5)
    ),
    #`Concentration (g/L)` = `Amount (g)` * Ratio %>% round(1),
    `Stock conc. (M)` = (`Amount (g)` / `Molar mass (g/mole)` * Ratio / `Water (L)`)  %>% round(5),
    `Final volume (mL)` = c(
        1, 1, 2, 1, 2.5,
        rep(0.1, 5)
    )
) %>%
    select(-Ratio)

ft <- tt %>%
    mutate(across(everything(), as.character)) %>%
    flextable() %>%
    colformat_md() %>%
    merge_at(i = 1:5, j = 1, part = "body") %>%
    merge_at(i = 6:10, j = 1, part = "body") %>%
    width(j = 2, width = 1.5) %>%
    width(j = 3, width = 1) %>%
    width(j = 4:5, width = 1) %>%
    width(j = 7, width = 1.5) %>%
    hline(i = 5, j = NULL, border = NULL, part = "body") %>%
    hline(i = 10, j = NULL, border = NULL, part = "body") %>%
    fix_border_issues() %>%
    # Add footnote
    footnote(i = 5, j = 2, ref_symbols = "a",
             value = as_paragraph("500 mL of Fe-EDTA is prepared by mixing two solutions at 50C: 1.4g of FeSO", as_sub("4")," in 0.25L of water and 1.85g of Na", as_sub("2"),"EDTA in 0.25L of water.")) %>%
    footnote(i = 5, j = 2, ref_symbols = "b",
             value = as_paragraph("Fe-EDTA is light sensitive. Make sure to wrap the bottle in tin foil.")) %>%
    font(fontname = "Time New Roman")

#0.02 M Fe-EDTA
# Stock     |  g   | water (L) | Concentration
# ----------|------|-----------|--------------
# FeSO~4~   | 1.4  |  0.25     | 0.1 M
# Na~2~EDTA | 1.85 |  0.25     | 0.1 M
#
# - [ ] mix at 50C to make 500 mL iron versenate.
# - [ ] Fe-EDTA is light sensitive, maje sure to wrap the bottle in tin foil.
#

save_as_image(ft, here::here("plots/TableS3.png"), webshot = "webshot2")




