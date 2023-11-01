cd
source ~/.zshrc


# Check
mamba activate scoary
mamba env list

scoary -g <gene_presence_absence.csv> -t <traits.csv>
