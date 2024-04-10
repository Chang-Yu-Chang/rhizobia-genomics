#' This script analyzes the trait data

renv::load()
library(tidyverse)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(broom.mixed) # for tidy up lme4
source(here::here("metadata.R"))

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
options(contrasts=c("contr.sum", "contr.poly"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_data, "phenotypes_analysis/symbiosis/plants.csv"))
gts <- read_csv(paste0(folder_data, "phenotypes_analysis/growth/gts.csv"))


#
