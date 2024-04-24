#' Plot the difference between sites in tmax, tmin

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
source(here::here("metadata.R"))

sites <- read_csv(paste0(folder_data, "phenotypes_analysis/sites/sites.csv"))
diff_vars <- read_csv(paste0(folder_data, "phenotypes_analysis/sites/diff_vars.csv"))

# Barplots for temperature contrast between populations
compute_mean <- function(diff_vars, pop) {
    diff_vars %>%
        filter(population == pop) %>%
        group_by(variable, yday) %>%
        summarize(diff_var = mean(diff_var))
}
plot_box <- function(diff_var_i) {
    diff_var_i %>%
        ggplot() +
        geom_hline(yintercept = 0, color = "maroon", alpha = 0.3) +
        geom_boxplot(aes(x = variable, y = diff_var), outlier.size = 0, outlier.color = NA, fill = "grey90") +
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

diff_var1 <- compute_mean(diff_vars, "VA")
diff_var2 <- compute_mean(diff_vars, "PA")

diff_vars %>%
    filter(population == "VA") %>%
    group_by(variable,resample) %>%
    summarize(diff_var = mean(diff_var)) %>%
    filter(variable == "tmax_deg_c") %>%
    pull(diff_var) %>%
    range()

set.seed(1)
p1 <- plot_box(diff_var1) +
    scale_y_continuous(limits = c(-2, 7), breaks = -2:7, expand = c(0,.1)) +
    ylab(expression(mean ~ "[" ~ t ~ "("~L~")" - t~ "("~H~")" ~ "]")) +
    ggtitle("Mountain sites")
p2 <- plot_box(diff_var2) +
    scale_y_continuous(limits = c(-2, 3), breaks = -2:7, expand = c(0,.1)) +
    ylab(expression(mean ~ "[" ~ t ~ "("~U~")" - t~ "("~S~")" ~ "]")) +
    ggtitle("City sites")

## STATS
t.test(diff_var1[diff_var1$variable == "tmax_deg_c",]$diff_var)
t.test(diff_var1[diff_var1$variable == "tmin_deg_c",]$diff_var)
t.test(diff_var2[diff_var2$variable == "tmax_deg_c",]$diff_var)
t.test(diff_var2[diff_var2$variable == "tmin_deg_c",]$diff_var)


#
p <- plot_grid(p1, p2, nrow = 1, align = "hv", axis = "lrbt", labels = c("A", "B"), scale = 0.9) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS2.png"), p, width = 6, height = 3)
