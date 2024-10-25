#' This script performs MK test on each gene

library(tidyverse)
library(iMKT)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

set_name = "elev_med"
#set_name = "urbn_mel"
tt <- read_gpas(set_name)

tt$list_sccg
# ff <- read_fsts(set_name)

# Check numbers
list_msas <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/msa_with_outgroup")) %>% str_subset(".fasta$")
length(list_msas)
list_msas <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/msa_with_outgroup_trimmed")) %>% str_subset(".fasta$")
length(list_msas)
#list_dafs <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/tables/")) %>% str_subset(".daf")
list_divs <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/tables/")) %>% str_subset(".div")
length(list_dafs)
length(list_divs)


# Filter for dafs and divs that are in the list of sccg
folder_tables <- paste0(folder_data, "genomics_analysis/mktest/", set_name, "/tables/")
tb <- tibble(div_file = list_divs) %>%
    mutate(gene = str_remove(div_file, ".div")) %>%
    mutate(daf_file = str_replace(div_file, ".div", ".daf")) %>%
    left_join(mutate(tt$list_sccg, is_sccg = T)) %>%
    filter(is_sccg)

tbs <- tb %>%
    mutate(
        daf_tb = map(daf_file, ~read_delim(paste0(folder_tables, .x), show_col_types = F)),
        div_tb = map(div_file, ~read_delim(paste0(folder_tables, .x), show_col_types = F))
    )


tbs_mktest <- tbs %>%
    # Some has 0 row in the daf? Does this mean lack of polymorphism?
    filter(map_lgl(daf_tb, ~nrow(.x)!=0)) %>%
    # Remove daf P0 all == 0
    filter(map_lgl(daf_tb, ~!all(.x$P0==0) & !all(.x$Pi==0))) %>%
    mutate(
        mktest = map2(daf_tb, div_tb, ~standardMKT(.x, .y)),
        mktest_alpha = map_dbl(mktest, ~`[[`(.x, 1)),
        mktest_p = map_dbl(mktest, ~`[[`(.x, 2)),
        mktest_mktable = map(mktest, ~`[[`(.x, 3)),
        mktest_divmetric = map(mktest, ~`[[`(.x, 4)),
    )

tbs_mktest %>%
    #unnest(mktest_divmetric)
    filter(mktest_p < 0.05) %>%
    select(gene, mktest_alpha) %>%
    view




gene = "zwf"
#gene = "accA"
gene = "aatB_2"
daf <- read_delim(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/tables/", gene, ".daf"))
div <- read_delim(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/tables/", gene, ".div"))
standardMKT(daf, div) %>% `[[`(1)

# myDafData
# myDivergenceData
# standardMKT(myDafData, myDivergenceData)

