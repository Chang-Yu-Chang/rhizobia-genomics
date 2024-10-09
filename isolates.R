
renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates %>%
    view

# sanger label id
read_csv(paste0(folder_data, "raw/rhizobia/02-sequencing/isolates_for_seq.csv")) %>%
    clean_names() %>%
    select(exp_id, sample_name) %>%
    view

# sanger taxonomy
read_csv(paste0(folder_data, "temp/deprecated/01-isolates_RDP.csv")) %>%
    clean_names() %>%
    select(exp_id, id, family, genus, family_score, genus_score) %>%
    view

# sanger rhizobia. list for wgs
read_csv(paste0(folder_data, "temp/deprecated/01-isolates_rhizo.csv")) %>%
    clean_names() %>%
    select(exp_id, id, family, genus, family_score, genus_score) %>%
    view

# wgs taxonomy
read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv")) %>%
    mutate(genome_id = factor(genome_id, paste0("g", 1:100))) %>%
    arrange(genome_id) %>%
    view


# VA lupulina
read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv")) %>%
    clean_names() %>%
    distinct(rhizobia) %>%
    view

# VA sativa
read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv")) %>%
    clean_names() %>%
    distinct(rhizobia) %>%
    view


# PA lupulina
read_csv(paste0(folder_data, "raw/SymbiosisInSoilData_S24.csv")) %>%
    clean_names() %>%
    distinct(rhizobia_strain) %>%
    view

# PA sativa
read_csv(paste0(folder_data, "raw/BIOL1102_PooledData.csv")) %>%
    clean_names() %>%
    distinct(rhizobia_strain, location) %>%
    view


# All
tb <- read_csv(paste0(folder_data, "mapping/isolates_all.csv"))
tb %>%
    filter(exp_id != 50) %>%
    filter(!is.na(exp_lup) | !is.na(exp_sativa)) %>%
    group_by(site_group) %>%
    count()

table(tb$exp_lup, tb$exp_sativa)



