#'

library(tidyverse)
library(cowplot)
#library(seqinr)
source(here::here("analysis/00-metadata.R"))

ant <- read_delim(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_1/05-bakta/consensus.tsv"), skip = 2)


dim(ant)

ant %>%
    drop_na(Gene) %>%
    view
