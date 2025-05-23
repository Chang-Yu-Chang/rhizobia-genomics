#' This script generates the table of ncbi accessions

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Table
en <- read_csv(paste0(folder_data, "raw/ensifer_ncbi.csv"), col_names = F)
en <- en[,1:3]
colnames(en) <- c("Accession", "Species", "Strain")
en <- en %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())


plasmids <- tibble(
    species = c(rep("meliloti", 3*4+1), rep("medicae", 3*4+1), "canadensis", "adhaerens"),
    strain = c(rep(c("usda1106", "em1021", "em1022"), each = 3), rep(c("usda1021", "wsm419", "wsm1115", "su277"), each = 4), "su277", "t173", "corn53"),
    replicon_type = c(rep(c("chromosome", "pSymA", "pSymB"), 3), rep(c("chromosome", "pSymA", "pSymB", "pAcce"), 4), "pAcce", "chromosome", "chromosome"),
    replicon = c(
        "chromosome", "psymA", "psymB",
        "chromosome", "pSymA", "pSymB",
        "chromosome", "pA", "pB",
        "chromosome", "psymA", "psymB", "accessoryA",
        "chromosome", "pSMED02", "pSMED01", "pSMED03",
        "chromosome", "pWSM1115_2", "pWSM1115_1", "pWSM1115_3",
        "chromosome", "pSU277_2", "pSU277_1", "pSU277_3", "pSU277_4",
        "chromosome", "chromosome"
    )
) %>%
    mutate(species = paste0("Sinorhizobium ", species)) %>%
    select(-species) %>%
    pivot_wider(id_cols = strain, names_from = replicon_type, values_from = replicon, values_fn = list) %>%
    mutate(
        chromosome = map_chr(chromosome, paste),
        pSymA = map_chr(pSymA, ~ if (is.null(.x)) "" else paste(.x, collapse = ", ")),
        pSymB = map_chr(pSymB, ~ if (is.null(.x)) "" else paste(.x, collapse = ", ")),
        pAcce = map_chr(pAcce, ~ if (is.null(.x)) "" else paste(.x, collapse = ", "))
    ) %>%
    rename(Strain = strain)


ft <- en %>%
    left_join(plasmids) %>%
    select(-chromosome) %>%
    flextable() %>%
    autofit() %>%
    merge_v(j = "Species") %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = c(1:4,8:11,13,17, 19)) %>%
    align(part = "all", align = "center") %>%
    valign(j = "Species", valign = "top") %>%
    style(j = "Species", pr_t = fp_text_default(italic = T)) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

#save_as_html(ft, path = here::here("plots/TabS1.html"))
save_as_image(ft, path = here::here("plots/TabS1.png"), res = 300)
