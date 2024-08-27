#' This script


# Cophenetic permutation test ----
source(here::here("forposter/cophenetic.R"))
tbp <- tb1 %>%
    mutate(id = factor(id, 1:12)) %>%
    unnest(distances_pm) %>%
    compute_percentiles()
plot_cophenetic <- function (tbp, tb_obs) {
    tbp %>%
        ggplot() +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, fill = alpha("gold", .2)) +
        geom_segment(aes(x = id, xend = id, y = p05, yend = p95), linewidth = 1,  arrow = arrow(length = unit(3, "mm"), angle = 90, ends = "both")) +
        geom_point(aes(x = id, y = p50, color = "95% CIs"), shape = 3, stroke = 1, size = 2) +
        geom_point(data = tb_obs, aes(x = id, y = distances_obs, color = "observation"), shape = 21, stroke = 2, size = 2) +
        scale_x_discrete(breaks = 1:12, labels = tb1$replicon_type) +
        scale_y_continuous(breaks = seq(0,1.5,0.5), limits = c(0, 1.5), minor_breaks = seq(0,1.5,0.1)) +
        scale_color_manual(values = c("observation" = "maroon", "95% CIs" = "black"), name = NULL) +
        facet_grid(.~ feature, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.8, 0.2),
            legend.background = element_rect(color = "black", fill = "white"),
            legend.box.background = element_rect(color = NA, fill = NA),
            legend.key = element_rect(color = NA, fill = NA),
            legend.spacing.y = unit(10,"mm"),
            legend.text = element_text(size = 8),
            legend.key.size = unit(5, "mm"),
            strip.background = element_blank(),
            strip.text = element_text(angle = 0, hjust = 0, vjust = 0, size = 8),
            strip.clip = "off",
            panel.spacing.x = unit(1, "mm"),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(size = 8, angle = 20, hjust = 1),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 15),
            plot.margin = unit(c(10,10,5,5), "mm"),
            plot.background = element_blank()
        ) +
        guides() +
        labs(x = "", y = expression(bar(d)["within"]/bar(d)["between"]))
}

p1 <- tbp %>%
    filter(population == "VA") %>%
    plot_cophenetic(filter(tb_obs, population == "VA")) +
    guides(color = "none") +
    labs(title = "Elevation gradient")
p2 <- tbp %>%
    filter(population == "PA") %>%
    plot_cophenetic(filter(tb_obs, population == "PA")) +
    labs(title = "Urbanization gradient")

p_coph <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb", label_y = .9, labels = c("C", "D"))

# Combine ----
p <- plot_grid(p_cong, p_coph, ncol = 1) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig5.png"), p, width = 10, height = 8)

