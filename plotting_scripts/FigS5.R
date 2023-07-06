#' Correlation between the foreground pixel number and root weight

library(tidyverse)
library(factoextra) # for plotting pca eclipse
source(here::here("analysis/00-metadata.R"))


treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
tt <- treatments %>%
    drop_na(c(root_weight_mg, network_area_px2))

# root biomass
p1 <- tt %>%
    mutate(network_area_px2 = network_area_px2/1e5) %>%
    ggplot() +
    geom_point(aes(x = root_weight_mg, y = network_area_px2), size = 2, shape = 21) +
    geom_smooth(aes(x = root_weight_mg, y = network_area_px2), method = "lm") +
    annotate("text", label = paste0("n=", nrow(tt)), x = -Inf, y = Inf, hjust = -1, vjust = 2) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "root weight (mg)", y = expression(paste("# of foreground pixels", (10^5))))


ggsave(here::here("plots/FigS5.png"), p, width = 4, height = 4)

cor.test(treatments$root_weight_mg, treatments$network_area_px2, method = "pearson") # cor = 0.901, p < 0.001
model <- lm(network_area_px2 ~ root_weight_mg, data = treatments)
summary(model)

