cd
source ~/.zshrc

# Check genome quality
mamba activate busco
mamba env list


medaka_consensus=$1
busco_folder=$2

busco -i $medaka_consensus --out_path $busco_folder -o "results" -m genome --auto-lineage-prok
# `-i` input sequence file in FASTA format
# `-out_path` output folder name
# `-o` output project name
# `-m genome` BUSCO analysis mode. We use genome mode
# `--auto-lineage-prok` Automated lineage selection


# busco -i ~/MGJW/problem_set3/out_S56/contigs.fa -l bacteria_odb10 -o busco_op -m genome
# busco -i ~/MGJW/problem_set1/fasta/genome4.fasta -l bacteria_odb10 -o busco_op2 -m genome
# busco -i ~/MGJW/problem_set1/fasta/genome4.fasta -l clostridiales_odb10 -o busco_op3 -m genome
# busco -i ~/MGJW/problem_set3/out_S56/contigs.fa -l pseudomonadales_odb10 -o busco_op4 -m genome
# Check genome completeness and contamination
# mamba activate checkm
# mamba env list
#
# checkm_results=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/06-checkm/checkm-results.tsv
# bin_input=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka
# checkm_dir=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/06-checkm
#
# checkm lineage_wf --tab_table --file "$checkm_results" --thread 12 --pplacer_threads 12 --extension fasta "$bin_input" "$checkm_dir"
# `lineage_wf` runs tree, lineage_set, analyze, qa
# `--tab_table` for tabular outputs, print a tab-separated values table instead of a table formatted for console output
# `--file` print output to file (default: stdout)
# `--thread 12` number of threads (default: 1)
# `--pplacer_threads 12` number of threads used by pplacer (memory usage increases linearly with additional threads) (default: 1)
# `--extension fasta`  extension of bins (other files in directory are ignored) (default: fna)
