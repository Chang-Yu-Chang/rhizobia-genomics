#' This script is to have a overview on the read statistics

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))


list_g <- paste0(rep("Chang_Q5C_results", 20), "/Chang_Q5C_", 1:20, "/")
list_g[11] <- "Chang_Q5C_results_repeated/Chang_Q5C_11/"
list_g[18] <- "Chang_Q5C_results_repeated/Chang_Q5C_18/"


for (i in 1:20) {
    raw_phred <- read_table(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "01-filtlong/raw_phred.txt"), col_names = c("name", "phred", "length"))
    p <- raw_phred %>%
        mutate(phred = phred - 33) %>%
        ggplot() +
        geom_point(aes(x = length, y = phred), size = 0.2, alpha = 0.2) +
        theme_classic() +
        theme(
            panel.grid.major = element_line(color = "grey90", linewidth = 0.1)
        ) +
        guides() +
        labs()

    ggsave(paste0(folder_data, "temp/31-01-raw_reads_g", sprintf("%02d", i),".png"), plot = p, width = 5, height = 5)

}

#filtered_phred <- read_table(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/01-filtered_phred.txt"), col_names = c("name", "phred", "length"))
#downsampled2_phred <- read_table(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye/03-downsampled2_phred.txt"), col_names = c("name", "phred", "length"))






























raw_phred %>%
    arrange(desc(length)) %>%
    mutate(phred = phred - 33) %>%
    summarize(median_phred = median(phred))

raw_phred %>%
    mutate(phred = phred - 33) %>%
    arrange(length)

median(raw_phred$length)
range(raw_phred$phred)
range(raw_phred$length)


raw_phred %>%
    arrange(desc(length)) %>%
    mutate(phred = phred - 33) %>%
    ggplot() +
    geom_point(aes(x = length, y = phred), size = 0.2, alpha = 0.2) +
    theme_classic() +
    theme() +
    guides() +
    labs()



filtered_phred  %>%
    arrange(desc(length)) %>%
    mutate(phred = phred - 33) %>%
    ggplot() +
    geom_point(aes(x = length, y = phred), size = 0.2, alpha = 0.2) +
    theme_classic() +
    theme() +
    guides() +
    labs()




test_read <- raw_phred[1,]
nchar(test_read$phred)
nchar(test_read$length)
p_cov <- function (qq) 10^(-qq/10)

p_cov(30)

median(raw_phred$phred) -33
median(filtered_phred$phred) -33
median(downsampled2_phred$phred) - 33

tt <- raw_phred %>%
    mutate(length)
    arrange(desc(length))
x <- tt$phred[1]

str_sub(x, 1, 10) %>% charToRaw() %>% as.numeric
charToRaw(x) %>% as.numeric() %>% mean()

str_detect(x, "J")


test_read <- tibble(asc = strsplit(x, "")[[1]]) %>%
    rowwise() %>%
    mutate(phred = as.numeric(charToRaw(asc)))

mean(test_read$phred)

test_read %>%
    mutate(phred = phred - 64) %>%
    ggplot() +
    geom_histogram(aes(x = phred)) +
    theme_classic() +
    theme() +
    guides() +
    labs()


charToRaw("|") %>% as.numeric()
str_detect(x, "J")

x <- "feeffdbefc`\\KKX]_BBBB"
charToRaw(x)
as.numeric(charToRaw(x))


raw_phred %>%
    ggplot() +
    geom_histogram(aes(x = length)) +
    theme_classic() +
    theme() +
    guides() +
    labs()



filtered_phred %>%
    slice(1:10000) %>%
    mutate(phred = phred-33)
    ggplot() +
    geom_point(aes(x = length, y = phred), size = 0.2, alpha = 0.1) +
    theme_classic() +
    theme() +
    guides() +
    labs()
