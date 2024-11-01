# Plot growth rate pcs

# 2. Plot the PCAs of 30C ----
gts1 <- gtwl %>%
    filter(gradient == "elevation", temperature == "30c") %>%
    pivot_wider(names_from = trait, values_from = value) %>%
    select(population, genome_id, r, lag, maxOD)
gts2 <- gtwl %>%
    filter(gradient == "urbanization", temperature == "30c") %>%
    pivot_wider(names_from = trait, values_from = value) %>%
    select(population, genome_id, r, lag, maxOD)

pca_results <- list(elevation = prcomp(gts1[,-c(1,2)], scale. = TRUE),
                    urbanization = prcomp(gts2[,-c(1,2)], scale. = TRUE))
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

pcs <- bind_rows(
    as_tibble(pca_results$elevation$x) %>% mutate(population = gts1$population, genome_id = gts1$genome_id),
    as_tibble(pca_results$urbanization$x) %>% mutate(population = gts2$population, genome_id = gts2$genome_id)
) %>%
    left_join(distinct(isolates, gradient, population)) %>%
    mutate(gradient = factor(gradient, c("elevation", "urbanization")))

plot_pca <- function (pcs, grad) {
    pcsi <- pcs %>% filter(gradient == grad)
    dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
    mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = genome_id))
    pcsi %>%
        ggplot() +
        geom_point(aes(x = PC1, y = PC2, color = population), shape = 21, stroke = 1, size = 2) +
        stat_ellipse(aes(x = PC1, y = PC2, fill = population), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 5))) +
        geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
        geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
        scale_color_manual(values = population_colors) +
        scale_fill_manual(values = population_colors, name = "population") +
        scale_x_continuous(breaks = seq(-8,8,2)) +
        scale_y_continuous(breaks = seq(-8,8,2)) +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey10", fill = NA),
            legend.background = element_blank(),
            legend.position = "top",
            legend.key = element_blank(),
            legend.title = element_blank(),
            aspect.ratio = 1
        ) +
        guides(color = "none") +
        labs(x = paste0("PC1 (", get_pcvar(pca_results[grad][[1]])[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results[grad][[1]])[2], "%)"))

}
p_pca1 <- plot_pca(pcs, "elevation")
p_pca2 <- plot_pca(pcs, "urbanization")

p <- plot_grid(p_pca1, p_pca2, nrow = 1, axis = "lr", align = "v") + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "phenotypes/growth/02-30c_pca.png"), p, width = 8, height = 4)

# Stats
pcsi <- pcs %>% filter(gradient == "elevation")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = genome_id))

pcsi <- pcs %>% filter(gradient == "urbanization")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = genome_id))

