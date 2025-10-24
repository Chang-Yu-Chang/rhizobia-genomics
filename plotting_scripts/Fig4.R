#' Plot IBD and compute Mantel tests for core (SCCG) and accessory (GCV) Dxy

library(tidyverse)
library(cowplot)
library(ggtext)  # for markdown in facet labels
source(here::here("metadata.R"))

# input ----
sccg_join <- read_csv(file.path(folder_genomics, "ibd/sccg_join.csv"), show_col_types = FALSE)
gcv_join  <- read_csv(file.path(folder_genomics, "ibd/gcv_join.csv"),  show_col_types = FALSE)
mantel_results <- read_csv(file.path(folder_genomics, "ibd/mantel_results.csv"), show_col_types = FALSE)

# relabel groups ----
label_map <- c(
    "meliloti_PA" = "*S.* *meliloti* in Pennsylvania",
    "medicae_VA"  = "*S.* *medicae* in Virginia"
)
color_map <- tibble(
    group = factor(label_map, label_map),
    color = species_colors[c(4,3)]
)

sccg_join <- sccg_join %>% mutate(group = recode(group, !!!label_map) %>% factor(label_map))
gcv_join  <- gcv_join  %>% mutate(group = recode(group, !!!label_map) %>% factor(label_map))
mantel_results <- mantel_results %>% mutate(group = recode(group, !!!label_map)%>% factor(label_map))

# join ----
mantel_results <- mantel_results %>%
    left_join(
        bind_rows(
            sccg_join %>%
                distinct(replicon, group, total_sites = total_sites_poly) %>%
                mutate(type = "Core genome (sequence similarity)"),
            gcv_join %>%
                distinct(replicon, group, total_sites = total_genes) %>%
                mutate(type = "Accessory genome (gene content similarity)")
        )
    ) %>%
    mutate(ast = map(p, detect_sig))

# plotting helper ----
plot_ibd_panel <- function(df, ylab, site_name, ty) {
    ggplot(df) +
        geom_rect(data = color_map, aes(fill = group), alpha = .2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        geom_point(aes(x = dist_geo_km, y = Dxy), shape = 21, color = "black", fill = NA, size = 1.5) +
        geom_smooth(aes(x = dist_geo_km, y = Dxy), method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
        geom_text(
            data = filter(mantel_results, str_detect(type, ty)),
            aes(label = paste0(total_sites, " ", site_name, "s, r²=", round(r, 3), " ", ast)),
            x = -Inf, y = Inf, hjust = -.05, vjust = 1.5, size = 2.5, color = "black"
        ) +
        facet_grid(group ~ replicon) +
        labs(y = ylab, x = "Geographic distance (km)") +
        scale_fill_manual(values = setNames(color_map$color, color_map$group)) +
        scale_y_continuous(expand = c(0, .1)) +
        theme_bw(base_size = 11) +
        theme(
            strip.background = element_rect(fill = NA, color = NA),
            strip.text.x = element_text(size = 10, face = "bold"),
            strip.text.y = element_markdown(face = "plain", size = 8),
            panel.grid = element_blank(),
            # axis.title.y = element_text(face = "bold"),
            plot.margin = margin(5, 5, 5, 5),
            legend.position = "none"
        )
}

# plot ----
p_core <- plot_ibd_panel(sccg_join, "Core genome Dxy", "SNP", "Core")
p_gcv  <- plot_ibd_panel(gcv_join, "Accessory genome Dxy", "gene", "Accessory")

# combine panels ----
p <- plot_grid(p_core, p_gcv, ncol = 1, scale = .95, labels = c("A", "B")) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig4.png"), p, width = 6, height = 8, dpi = 300)

