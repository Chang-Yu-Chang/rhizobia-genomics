cd
source ~/.zshrc

# Check
mamba activate roary
mamba env list

folder_temp_result=/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus
cd $folder_temp_result
mkdir -p summary/42-roary

# Copy the gff files to a folder
cd summary/42-roary
mkdir -p gff

for g in {1..19}
do
    echo g$g
    cp $folder_temp_result/Chang_Q5C_$g/10-prokka/annotated.gff $folder_temp_result/summary/42-roary/gff/annotated_g$g.gff
done

# Run Roary
cd $folder_temp_result/summary/42-roary
roary -f roary1 -e -n -i 95 -cd 99 -s gff/*.gff
# `-f STR` output directory
# `-e` create a multiFASTA alignment of core genes using PRANK
# `-n` fast core gene alignment with MAFFT
# `-i 95` minimum percentage identity for blastq [95]
# `-cd FLOAT` percentage of isolates in a gene must be in to be core [99]
# `-r` create R plots
# `-s` dont split paralogs
# `*.gff` input

# Run Roary again, changing the cut-off scores (-i) to 75, then again with 50. How do the results change?
cd $folder_temp_result/summary/42-roary
roary -f roary2 -e -n -i 75 -cd 99 -s gff/*.gff

cd $folder_temp_result/summary/42-roary
roary -f roary3 -e -n -i 50 -cd 99 -s gff/*.gff

cd $folder_temp_result/summary/42-roary
roary -f roary4 -e -n -i 95 -cd 90 -s gff/*.gff


# Run a fast tree
FastTree -nt -gtr core_gene_alignment.aln > my_tree.newick
# QUick plot roary results. Make sure to specify the correct directory of roary_plots.py
python3 roary_plots.py my_tree.newick gene_presence_absence.csv

"finish the description of pangenome"

# Check the results
# for rr in roary1 roary2 roary3 roary4
# do
#     echo $rr
#     cat $rr/summary_statistics.txt
# done



cd roary1
# Union of genes found in isolates. The output is a table
query_pan_genome -a union ../gff/*.gff
# Intersection of genes found in isolates (core genes)
query_pan_genome -a intersection ../gff/*.gff
# Complement of genes found in isolates (accessory genes)
query_pan_genome -a complement ../gff/*.gff
