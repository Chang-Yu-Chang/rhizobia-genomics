#' This script performs MK test on each gene

library(tidyverse)
library(iMKT)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

check_files <- function (set_name, ref) {
    # Check numbers
    list_msas <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/msa_with_outgroup_trimmed")) %>% str_subset(".fasta$")
    cat("\nTrimmed msa with outgroup: ",length(list_msas))
    list_dafs <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/tables/")) %>% str_subset(".daf")
    cat("\nDAF: ", length(list_dafs))
    list_divs <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/tables/")) %>% str_subset(".div")
    cat("\nDIV: ", length(list_divs))
}
aggregate_tbs <- function (set_name, ref) {
    tt <- read_gpas(set_name)
    list_divs <- list.files(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/tables/")) %>% str_subset(".div")

    # Filter for dafs and divs that are in the list of sccg
    folder_tables <- paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/tables/")
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

    cat("\nNumber of genes in MK test: ", nrow(tbs_mktest))

    mktest <- tbs_mktest %>% select(gene, mktest_alpha, mktest_p)
    write_csv(mktest, paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/mktests.csv"))
    #return(mktest)
}

#check_files("elev_med", "ngr234")
aggregate_tbs("elev_med", "ngr234")
aggregate_tbs("urbn_mel", "ngr234")
#aggregate_tbs("elev_med", "em1021")
#aggregate_tbs("urbn_mel", "wsm419")
