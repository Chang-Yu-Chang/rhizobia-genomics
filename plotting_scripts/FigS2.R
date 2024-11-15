#' Plot the difference between sites in tmax, tmin

library(tidyverse)
library(cowplot)
source(here::here("metadata.R"))

sites <- read_csv(paste0(folder_data, "phenotypes/sites/sites.csv"))
diff_vars <- read_csv(paste0(folder_data, "phenotypes/sites/diff_vars.csv"))

# Barplots for temperature contrast between populations
compute_mean <- function(diff_vars, gra) {
    diff_vars %>%
        filter(gradient == gra) %>%
        group_by(variable, yday) %>%
        summarize(diff_var = mean(diff_var))
}
plot_box <- function(diff_var_i) {
    diff_var_i %>%
        ggplot() +
        geom_hline(yintercept = 0, color = "maroon", alpha = 0.3) +
        geom_boxplot(aes(x = variable, y = diff_var), outlier.size = 0, outlier.color = NA, fill = "grey90") +
        geom_jitter(aes(x = variable, y = diff_var), shape = 21, size = 2, color = "black", width = 0.2, alpha = 0.5) +
        scale_x_discrete(breaks = c("tmax_deg_c", "tmin_deg_c"), labels = c(expression(t[max]), expression(t[min]))) +
        theme_classic() +
        theme(
            panel.grid.major.y = element_line(color = "grey90", linewidth = .5, linetype = 1),
            panel.grid.minor.y = element_line(color = "grey90", linewidth = .5, linetype = 2),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title.y = element_text(size = 10)
        ) +
        guides() +
        labs()
}

diff_var1 <- compute_mean(diff_vars, "elevation")
diff_var2 <- compute_mean(diff_vars, "urbanization")

diff_vars %>%
    filter(gradient == "elevation") %>%
    group_by(variable,resample) %>%
    summarize(diff_var = mean(diff_var)) %>%
    filter(variable == "tmax_deg_c") %>%
    pull(diff_var) %>%
    range()

set.seed(1)
p1 <- plot_box(diff_var1) +
    scale_y_continuous(limits = c(-2, 7), breaks = -2:7, expand = c(0,.1)) +
    ylab(expression(mean ~ "[" ~ t ~ "("~L~")" - t~ "("~H~")" ~ "] (oC)")) +
    ggtitle("Elevation gradient")
p2 <- plot_box(diff_var2) +
    scale_y_continuous(limits = c(-2, 3), breaks = -2:7, expand = c(0,.1)) +
    ylab(expression(mean ~ "[" ~ t ~ "("~U~")" - t~ "("~S~")" ~ "] (oC)")) +
    ggtitle("Urbanization gradient")

## STATS
t.test(diff_var1[diff_var1$variable == "tmax_deg_c",]$diff_var)
t.test(diff_var1[diff_var1$variable == "tmin_deg_c",]$diff_var)
t.test(diff_var2[diff_var2$variable == "tmax_deg_c",]$diff_var)
t.test(diff_var2[diff_var2$variable == "tmin_deg_c",]$diff_var)


#
p <- plot_grid(p1, p2, nrow = 1, align = "hv", axis = "lrbt", labels = c("A", "B"), scale = 0.9) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS2.png"), p, width = 5, height = 4)



