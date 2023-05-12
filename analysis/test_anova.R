#' Practice ANOVA

library(tidyverse)

data(PlantGrowth)
PlantGrowth <- as_tibble(PlantGrowth)
levels(PlantGrowth$group)

PlantGrowth %>%
    ggplot() +
    geom_point(aes(x = group, y = weight)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

#
mod <- aov(weight ~ group, data = PlantGrowth)
summary(mod)

# Sum of squares total
ss_total <- sum((PlantGrowth$weight - mean(PlantGrowth$weight))^2) # 14.25843

# Sum of squares explained. Between group
ss_explained <- sum((PlantGrowth$group_mean - mean(PlantGrowth$weight))^2) # 3.76634

# Sum of squares residual. Within group
group_means <- tapply(PlantGrowth$weight, PlantGrowth$group, mean)
group_means <- tibble(group = names(group_means), group_mean = group_means)
PlantGrowth <- PlantGrowth %>% left_join(group_means)
ss_residual <- sum((PlantGrowth$weight - PlantGrowth$group_mean)^2) # 10.492

# df_explained = df of the explained part = number of group = 1
n_groups <- length(unique(PlantGrowth$group))
df_explained = n_groups - 1 # 2
# df_residual = df of the residual = number of observations - number of groups
n_obs <- nrow(PlantGrowth)
df_residual <- n_obs - n_groups # 27

# Mean squares explained
ms_explained <- sum((PlantGrowth$group_mean - mean(PlantGrowth$weight))^2) / 2 # 1.8832

# Mean squares residual
ms_residual <- sum((PlantGrowth$weight - PlantGrowth$group_mean)^2) / 27 # 0.3886

# F value = MS_explained / residual
f <- ms_explained / ms_residual # 4.846

# p value
p_value <- 1 - pf(f, df_explained, df_residual) # 0.0159


1 - pf(1, df_explained, df_residual)








