cd
source ~/.zshrc

cd ~/Desktop/lab/local-adaptation/shell_scripts
source 00-env_vars.sh
batch_name=$1
#batch_name="Chang_Q5C_results"
sample_id=$2
#sample_id="Chang_Q5C_2"
echo $sample_id
echo $folder_genomes

# Make folders and specify log
for wf in 01-reads_qc 02-guided_assembly 03-denovo_assembly 04-taxonomy_assignment 05-gene_annotation
do
mkdir -p "$folder_genomes/$sample_id/$wf"
done

# 1. Quality control raw reads
# Filter the worst 5% reads via filtlong
echo "01-filtlong"
zsh 01-filter_reads.sh \
    "$folder_raw/$batch_name/$sample_id/reads/raw_reads.fastq.gz" \
    "$folder_genomes/$sample_id/01-reads_qc/filtered_reads.fastq.gz"

# Extract the filter
zsh 01a-extract_reads.sh \
    "$folder_genomes/$sample_id/01-reads_qc/filtered_reads.fastq.gz" \
    "$folder_genomes/$sample_id/01-reads_qc/filtered_reads.txt"

# Plot the read data
Rscript 01b-plot_reads.R \
    "$folder_genomes/$sample_id/01-reads_qc/filtered_reads.txt"\
    "$folder_genomes/$sample_id/01-reads_qc/filtered_reads_qc.png"


# 2.

# # 0. nanoplot
# echo "0-nanoplot"
# raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
# nanoplot_folder="$folder_genomes/$sample_id/00-nanoplot"
#
# zsh 30-nanoplot.sh $raw_reads $nanoplot_folder &> $log00
#
# # 2. Draft genome via miniasm
# echo "2-miniasm"
# downsample_reads="$folder_genomes/$sample_id/02-miniasm/02-downsampled_reads.fastq"
# assemble_paf="$folder_genomes/$sample_id/02-miniasm/02-miniasm_pilot.paf.gz"
# assemble_gfa="$folder_genomes/$sample_id/02-miniasm/02-miniasm_pilot.gfa"
# assemble_fq="$folder_genomes/$sample_id//02-miniasm/02-miniasm_pilot.fasta"
#
# zsh 32-miniasm.sh $filtered_reads $downsample_reads $assemble_paf $assemble_gfa $assemble_fq &> $log02
#
# # 3. Assemble genome via flye
# echo "3-flye"
# filtered_reads="$folder_genomes/$sample_id/01-filtlong/01-filtered_reads.fastq.gz"
# assemble_fq="$folder_genomes/$sample_id/02-miniasm/02-miniasm_pilot.fasta"
# downsampled_reads2="$folder_genomes/$sample_id/03-flye/03-downsampled2_reads.fastq"
# genome_size=$(grep -v ">" "$assemble_fq" | grep -E -o "G|C|T|A|N" | wc -l)
# flye_folder="$folder_genomes/$sample_id/03-flye/"
#
# zsh 33-flye.sh $filtered_reads $assemble_fq $downsampled_reads2 $genome_size $flye_folder &> $log03
#
# # 4. Polish genome via medaka
# echo "4-medaka"
# raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
# draft_assembly="$folder_genomes/$sample_id/03-flye/assembly.fasta"
# medaka_folder="$folder_genomes/$sample_id/04-medaka"
#
# zsh 34-medaka.sh $raw_reads $draft_assembly $medaka_folder &> $log04
#
# # 5. Annotate genome via bakta
# echo "5-bakta"
# bakta_database="/Users/cychang/bioinformatics/bakta/db" # This database is mandatory and must be downloaded before annotation
# medaka_consensus="$folder_genomes/$sample_id/04-medaka/consensus.fasta"
# bakta_folder="$folder_genomes/$sample_id/05-bakta"
#
# zsh 35-bakta.sh $bakta_database $medaka_consensus $bakta_folder &> $log05
#
# # 6. Check genome quality via quast
# echo "6-quast"
# medaka_consensus="$folder_genomes/$sample_id/04-medaka/consensus.fasta"
# quast_folder="$folder_genomes/$sample_id/06-quast"
#
# zsh 36-quast.sh $medaka_consensus $quast_folder &> $log06
#
# # 7. Check genome quality via busco
# echo "7-busco"
# medaka_consensus="$folder_genomes/$sample_id/04-medaka/consensus.fasta"
# busco_folder="$folder_genomes/$sample_id/07-busco"
#
# zsh 37-busco.sh $medaka_consensus $busco_folder &> $log07
#
# # 8. Identify taxonomy via mash
# echo "8-mash"
# medaka_consensus="$folder_genomes/$sample_id/04-medaka/consensus.fasta"
# mash_folder="$folder_genomes/$sample_id/08-mash"
#
# zsh 38-mash.sh $medaka_consensus $mash_folder &> $log08
#
# # 9. Identify taxonomy via sourmash
# echo "9-sourmash"
# medaka_consensus="$folder_genomes/$sample_id/04-medaka/consensus.fasta"
# sourmash_folder="$folder_genomes/$sample_id/09-sourmash"
#
# zsh 39-sourmash.sh $medaka_consensus $sourmash_folder &> $log09
#
# # 10. prokka
# medaka_consensus="$folder_genomes/$sample_id/04-medaka/consensus.fasta"
# prokka_folder="$folder_genomes/$sample_id/10-prokka"
#
# zsh 40-prokka.sh $medaka_consensus $prokka_folder &> $log10

# 10. plasforest
# medaka_consensus="$folder_genomes/$sample_id/04-medaka/consensus.fasta"
# plasforest_folder="$folder_genomes/$sample_id/45-plasforest"
#
# zsh 45-plasforest.sh $medaka_consensus $plasforest_folder &> $log45



# # 1111. multiqc
# mamba activate multiqc
# cd folder_temp_result="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/$1/$2/"
# multiqc .
#
# # Remove all temporary files
#
#
# #
# #
# # Make folders
# mkdir -p "$folder_genomes/$sample_id"
# mkdir -p "$folder_genomes/$sample_id/logs"
# mkdir -p "$folder_genomes/$sample_id/00-nanoplot"
# mkdir -p "$folder_genomes/$sample_id/01-filtlong"
# mkdir -p "$folder_genomes/$sample_id/02-miniasm"
# mkdir -p "$folder_genomes/$sample_id/03-flye"
# mkdir -p "$folder_genomes/$sample_id/04-medaka"
# mkdir -p "$folder_genomes/$sample_id/05-bakta"
# mkdir -p "$folder_genomes/$sample_id/06-quast"
# mkdir -p "$folder_genomes/$sample_id/07-busco"
# mkdir -p "$folder_genomes/$sample_id/08-mash"
# mkdir -p "$folder_genomes/$sample_id/09-sourmash"
# mkdir -p "$folder_genomes/$sample_id/10-prokka"
# mkdir -p "$folder_genomes/$sample_id/45-plasforest"
#
# # Specify log directory
# log00="$folder_genomes/$sample_id/logs/00-nanoplot.log"
# log01="$folder_genomes/$sample_id/logs/01-filter.log"
# log02="$folder_genomes/$sample_id/logs/02-miniasm.log"
# log03="$folder_genomes/$sample_id/logs/03-flye.log"
# log04="$folder_genomes/$sample_id/logs/04-medaka.log"
# log05="$folder_genomes/$sample_id/logs/05-bakta.log"
# log06="$folder_genomes/$sample_id/logs/06-quast.log"
# log07="$folder_genomes/$sample_id/logs/07-busco.log"
# log08="$folder_genomes/$sample_id/logs/08-mash.log"
# log09="$folder_genomes/$sample_id/logs/09-sourmash.log"
# log10="$folder_genomes/$sample_id/logs/10-prokka.log"
# log45="$folder_genomes/$sample_id/logs/45-plasforest.log"
