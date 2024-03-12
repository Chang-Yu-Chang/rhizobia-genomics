#' This script generates figure 1

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(ggmap)
#library(ggsn) # add scale to map
library(ggsci)
library(RColorBrewer)
source(here::here("analysis/00-metadata.R"))

sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"), show_col_types = F)
diff_vars <- read_csv(paste0(folder_data, "temp/22-diff_vars.csv"), show_col_types = F)
population_colors <- c("PA" = "gold2", "VA" = "olivedrab")
site_group_colors <- c(`high elevation` = "#0C6291", `low elevation` = "#BF4342", `suburban` = "#0cc45f", `urban` = "#a642bf", control = "grey")


# Panel A. Cartoon for characterizing divergence in two populations
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png"))

# Panel B. Plot both populations in one map 
google_api <- "AIzaSyCmRu9YUWIQzq5U48yDeFPAiUYSft86pWw"
register_google(google_api)

# Sites
mlbs_center <- c(-80.55, 37.35)
mlbs_length <- 350
map_mlbs <- get_googlemap(center = mlbs_center, maptype = "terrain", size = c(mlbs_length, mlbs_length), style = "element:labels|visibility:off")
phila_center <- c(-75.25, 39.9)
phila_length <- 350
map_phila <- get_googlemap(center = phila_center, maptype = "terrain", size = c(phila_length, phila_length), style = "element:labels|visibility:off")

# VA map
dis <- 0.125
p2_1 <- ggmap(map_mlbs) +
    geom_point(data = filter(sites, population == "VA"), aes(x = longitude_dec, y = latitude_dec, color = site_group), shape = 21, size = 3, stroke = 1, fill = NA) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(limits = mlbs_center[1] + c(-dis, dis), breaks = ) +
    scale_y_continuous(limits = mlbs_center[2] + c(-dis, dis)) +
    coord_fixed() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.3, 0.85),
        legend.margin = margin(-5,1,1,1),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs(x = expression(longitude(~degree~W)), y = expression(latitude(~degree~N)))

## PA map
p2_2 <- ggmap(map_phila) +
    geom_point(data = filter(sites, population == "PA"), aes(x = longitude_dec, y = latitude_dec, color = site_group), shape = 21, size = 3, stroke = 1, fill = NA) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(limits = phila_center[1] + c(-dis, dis), breaks = phila_center[1] + c(0.1, 0, -0.1)) +
    scale_y_continuous(limits = phila_center[2] + c(-dis, dis), breaks = phila_center[2] + c(0.05, -0.05)) +
    coord_fixed() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.15),
        legend.key.size = unit(5, "mm"),
        legend.margin = margin(-5,1,1,1),
        legend.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs(x = expression(longitude(~degree~W)), y = expression(latitude(~degree~N)))

p2 <- plot_grid(p2_1, p2_2, ncol = 1, scale = 1, align = "vh", axis = "lrtb")


# Panel C. Barplots for temperature contrast between populations
#p3 <- ggdraw() + draw_text("placeholder")
compute_mean <- function(diff_vars, pop) {
    diff_vars %>%
        filter(population == pop) %>%
        group_by(variable, yday) %>%
        summarize(diff_var = mean(diff_var))
}

diff_var1 <- compute_mean(diff_vars, "VA") 
diff_var2 <- compute_mean(diff_vars, "PA") 

diff_vars %>%
    filter(population == "VA") %>%
    group_by(variable,resample) %>%
    summarize(diff_var = mean(diff_var)) %>%
    filter(variable == "tmax_deg_c") %>%
    pull(diff_var) %>%
    range()

plot_box <- function(diff_var_i) {
    diff_var_i %>%
        ggplot() +
        geom_hline(yintercept = 0, color = "black") +
        geom_boxplot(aes(x = variable, y = diff_var), outlier.size = 0) +
        geom_jitter(aes(x = variable, y = diff_var), shape = 21, size = 1, color = "black", width = 0.2, alpha = 0.3) +
        scale_x_discrete(breaks = c("tmax_deg_c", "tmin_deg_c"), labels = c(expression(t[max]), expression(t[min]))) +
        theme_classic() +
        theme(
            panel.grid.major.y = element_line(color = "grey90", linewidth = .5, linetype = 1),
            panel.grid.minor.y = element_line(color = "grey90", linewidth = .5, linetype = 2),
            axis.title.x = element_blank()
        ) +
        guides() +
        labs()
}

p3_1 <- plot_box(diff_var1) + 
    scale_y_continuous(limits = c(-2, 7), breaks = -2:7, expand = c(0,0)) + 
    ylab(expression(mean ~ "[" ~ t ~ "("~L~")" - t~ "("~H~")" ~ "]"))
p3_2 <- plot_box(diff_var2) +
    scale_y_continuous(limits = c(-2, 3), breaks = -2:7, expand = c(0,0)) +
    ylab(expression(mean ~ "[" ~ t ~ "("~U~")" - t~ "("~S~")" ~ "]"))
p3 <- plot_grid(p3_1, p3_2, nrow = 2, align = "hv", axis = "lrbt", scale = 0.9)

## STATS
t.test(diff_var1[diff_var1$variable == "tmax_deg_c",]$diff_var) 
t.test(diff_var1[diff_var1$variable == "tmin_deg_c",]$diff_var) 
t.test(diff_var2[diff_var2$variable == "tmax_deg_c",]$diff_var) 
t.test(diff_var2[diff_var2$variable == "tmin_deg_c",]$diff_var) 

# 
p <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1,.6,.5), labels = c("A", "B", "C")) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/Fig1.png"), p, width = 10, height = 5)
