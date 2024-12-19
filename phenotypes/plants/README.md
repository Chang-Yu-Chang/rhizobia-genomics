This folder contains scripts for handling symbiosis trait data

1. `clean_plants.R` cleans the variable names and binds the lupulina/sativa data tables into one `plants.csv`
2. `stat_trait_comparison.R` stats for pairwise population comparison
3. `plot_trait_comparison.R` plots the stats table
4. `stat_trait_all.R` stats for pairwise population comparison, using all sativa traits, including continuous and catagorical data

Scripts for computing effective size from symbiosis trait data

1. `compute_effectsize.R` computes effect size for each trait. There are three metrics: Cohen's d. Hedge's g, and partial eta squared
2. `plot_effectsize.R` plots the effect size diagram

Scripts for nitrogen treatment analysis

1. `stat_nitrogen_rn.R` performs permutation for reaction norm of nitrogen X population
2. `plot_nitrogen_rn.R` plots the stat table
