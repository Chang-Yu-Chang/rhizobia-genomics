#' This script plots the tradeoff figure

library(tidyverse)
library(janitor)
library(broom)
library(cowplot)
library(corrplot) # For comparing correlation
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>%
    filter(!genome_id %in% c("g20", "g28"))
gst_r <- read_csv(paste0(folder_data, "phenotypes/tradeoff/gst_r.csv"))
plants_r <- read_csv(paste0(folder_data, "phenotypes/tradeoff/plants_r.csv"))


# Plot the individual traits ----
plants_rs <- list()
plants_rs[[1]] <- filter(plants_r, exp_plant == "lupulina", exp_nitrogen == "without nitrogen")
plants_rs[[2]] <- filter(plants_r, exp_plant == "sativa", exp_nitrogen == "without nitrogen")
plants_rs[[3]] <- filter(plants_r, exp_plant == "sativa", exp_nitrogen == "with nitrogen")

list_treatments <- c("lupulina N-", "sativa N-", "sativa N+")

for (i in 1:3) {
    isolates_trait <- isolates %>%
        filter(!genome_id %in% c("g2", "g3", "g15")) %>%
        #filter(gradient == "urbanization") %>%
        left_join(plants_rs[[i]]) %>%
        left_join(gst_r) %>%
        drop_na(exp_nitrogen) %>%
        ungroup()

    m <- isolates_trait[,-c(1:8)]
    m <- m[,(colSums(is.na(m)) == 0 & colSums(m) != 0)] # Remove traits with NA
    mm <- cor(m)
    p_mat <- cor.mtest(m, conf.level = 0.95)
    tt <- min(which(str_detect(names(m), "^r_"))) # the index of first growth trait; for plotting

    png(paste0(folder_data, "phenotypes/tradeoff/0", i, ".png"), width = 20, height = 20, units = "cm", res = 300)

    corrplot(
        mm, diag = T, title = paste0(list_treatments[i], ", n_strains = ", nrow(isolates_trait)), mar = c(0,0,2,0),
        p.mat = p_mat$p, insig = "label_sig", sig.level = c(0.001, 0.01, 0.05), pch.cex = 1,
        type = "full", tl.col = "black", tl.cex = 1, cl.ratio = 0.1, tl.srt = 45, tl.offset = .5
    ) %>%
        corrRect(index = c(1, tt, ncol(m)))
    dev.off()
}


# Plot the PC1s ----
gst_pcs <- read_csv(paste0(folder_data, "phenotypes/tradeoff/gst_pcs.csv"))
gst_pc_imp <- read_csv(paste0(folder_data, "phenotypes/tradeoff/gst_pc_imp.csv"))
plants_pcs <- read_csv(paste0(folder_data, "phenotypes/tradeoff/plants_pcs.csv"))
plants_pc_imp <- read_csv(paste0(folder_data, "phenotypes/tradeoff/plants_pc_imp.csv"))

make_pc_lab <- function (pc, main_or_primary = 1) {
    paste0("PC", main_or_primary, " (", round(pc, 4)*100, "%)")
}
plot_pca <- function (pcs, imps) {
    pcs %>%
        ggplot() +
        geom_point(aes(x = PC1, y = PC2), shape = 21, size = 2) +
        #stat_ellipse(aes(x = PC1, y = PC2, fill = population), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
        geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
        geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
        #scale_fill_manual(values = population_colors, name = "population") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey10", fill = NA),
            plot.background = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            aspect.ratio = 1
        ) +
        guides(color = "none") +
        labs(x = make_pc_lab(imps$PC1, 1), y = make_pc_lab(imps$PC2, 2))
}

gst_pcs %>%
    left_join(isolates) %>%
    plot_pca(gst_pc_imp)


plants_pp <- list()
plants_pp[[1]] <- filter(plants_pcs, exp_plant == "lupulina", exp_nitrogen == "without nitrogen")
plants_pp[[2]] <- filter(plants_pcs, exp_plant == "sativa", exp_nitrogen == "without nitrogen")
plants_pp[[3]] <- filter(plants_pcs, exp_plant == "sativa", exp_nitrogen == "with nitrogen")

list_treatments <- c("lupulina N-", "sativa N-", "sativa N+")
isolates_pcs <- list()
for (i in 1:3) {
    isolates_trait <- isolates %>%
        filter(!genome_id %in% c("g2", "g3", "g15")) %>%
        #filter(gradient == "urbanization") %>%
        left_join(rename_with(plants_pp[[i]], .cols = starts_with("PC"), function (x) paste0("sym_", x))) %>%
        left_join(rename_with(gst_pcs, .cols = starts_with("PC"), function (x) paste0("gth_", x))) %>%
        left_join(gst_r) %>%
        left_join(plants_rs[[i]]) %>%
        drop_na(exp_nitrogen) %>%
        ungroup() %>%
        drop_na(sym_PC1, gth_PC1)

    isolates_pcs[[list_treatments[i]]] <- isolates_trait
}

plot_two_PC1s <- function (isolates_trait, x, y, use_PCs = T, by_gradient = T, sym_pc1 = 0, gst_pc1 = 0) {
    if (by_gradient == T) isolates_trait <- isolates_trait %>% mutate(gradient = factor(gradient, c("elevation", "urbanization")))

    p <- isolates_trait %>%
        ggplot() +
        geom_point(aes(x = {{x}}, y = {{y}}, color = population), shape = 21, size = 2, stroke = 1) +
        geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
        geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
        scale_color_manual(values = population_colors, name = "") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey10", fill = NA),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            aspect.ratio = 1,
            legend.position = "inside"
        ) +
        guides() +
        ggtitle(label = "", subtitle = paste0("n_strains = ", nrow(isolates_trait)))

    if (use_PCs == T) p <- p + labs(x = paste0("symbiosis ", make_pc_lab(sym_pc1, 1)), y = paste0("growth ", make_pc_lab(gst_pc1, 1)))
    if (by_gradient == T) p <- p + facet_grid(~gradient, drop = F)
    return(p)
}
make_cor_label <- function (cor_tidied) {
    paste0("Pearson's r=", round(cor_tidied$estimate, 2), ", p=", round(cor_tidied$p.value, 3))
}

# Plot by treatment ----
pc_results <- isolates_pcs %>%
    bind_rows(.id = "treatment") %>%
    #group_by(treatment, gradient) %>%
    group_by(treatment) %>%
    nest() %>%
    mutate(
        cor_result = map(data, ~cor.test(.x$sym_PC1, .x$gth_PC1)),
        cor_tidied = map(cor_result, ~tidy(.x))
    ) %>%
    unnest(cor_tidied) %>%
    ungroup %>%
    mutate(pc_imp_dummy_sym = c(1,2,3), pc_imp_dummy_gth = c(1,1,1))


p_list <- list()
for (i in 1:nrow(pc_results)) {
    p_list[[i]] <- plot_two_PC1s(
        pc_results$data[[i]], sym_PC1, maxOD_30c,
        use_PCs = T, plants_pc_imp$PC1[pc_results$pc_imp_dummy_sym[i]], gst_pc_imp$PC1[pc_results$pc_imp_dummy_gth[i]],
        by_gradient = F
    ) +
        guides(color = "none") +
        annotate("text", label = make_cor_label(pc_results[i,]), size = 3, x = Inf, y = -Inf, hjust = 1.05, vjust = -1.3) +
        ggtitle(pc_results$treatment[i])
}

p <- plot_grid(plotlist = p_list, ncol = 1, scale = 0.9, align = "v", axis = "rl") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_data, "phenotypes/tradeoff/pcs_treatment.png"), p, width = 3, height = 10)

# Plot by gradient and treatment ----
pc_results <- isolates_pcs %>%
    bind_rows(.id = "treatment") %>%
    group_by(treatment, gradient) %>%
    nest() %>%
    mutate(
        cor_result = map(data, ~cor.test(.x$sym_PC1, .x$gth_PC1)),
        cor_tidied = map(cor_result, ~tidy(.x))
    ) %>%
    unnest(cor_tidied) %>%
    ungroup %>%
    mutate(pc_imp_dummy_sym = c(1,1,2,2,3), pc_imp_dummy_gth = c(1,1,1,1,1))

p_list <- list()
for (i in 1:nrow(pc_results)) {
    p_list[[i]] <- plot_two_PC1s(
        pc_results$data[[i]], sym_PC1, maxOD_30c,
        use_PCs = T, plants_pc_imp$PC1[pc_results$pc_imp_dummy_sym[i]], gst_pc_imp$PC1[pc_results$pc_imp_dummy_gth[i]],
        by_gradient = F
    ) +
        guides(color = "none") +
        annotate("text", label = make_cor_label(pc_results[i,]), size = 3, x = Inf, y = -Inf, hjust = 1.05, vjust = -1.3) +
        ggtitle(pc_results$treatment[i])
}

p <- plot_grid(plotlist = p_list, ncol = 2, scale = 0.9, align = "v", axis = "rl") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_data, "phenotypes/tradeoff/pcs_gradient.png"), p, width = 6, height = 10)
#ggsave(paste0(folder_data, "phenotypes/tradeoff/pcs-0", i, ".png"), p, width = 6, height = 3)




# xx <- pc_results$data[[2]]
# tidy(cor.test(xx$sym_PC1, xx$gth_PC1))
# tidy(cor.test(xx$nodule_number, xx$maxOD_30c))
# xx$nodule_number
#
#
#
#
#
#
#
#
#
#
#
