cd
source ~/.zshrc


mamba activate scoary
mamba env list

gene_pa="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/42-roary/roary5/gene_presence_absence.csv"
traits="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/42-isolates.csv"

scoary -g $gene_pa -t $traits
"NEVER MIND. SCOARY only works for binary traits"
