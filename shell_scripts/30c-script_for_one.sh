cd
source ~/.zshrc

cd ~/Desktop/lab/local-adaptation/shell_scripts
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
folder_raw_result="/Users/cychang/Dropbox/lab/local-adaptation/data/raw/$1"
folder_temp_result="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/$1"

# Sample ID
sample_id=$2
mkdir -p "$folder_temp_result/$sample_id"
mkdir -p "$folder_temp_result/$sample_id/logs"

echo $folder_temp_result
echo $sample_id

# Make folders
mkdir -p "$folder_temp_result/$sample_id/01-filtlong"
mkdir -p "$folder_temp_result/$sample_id/02-miniasm"
mkdir -p "$folder_temp_result/$sample_id/03-flye"
mkdir -p "$folder_temp_result/$sample_id/04-medaka"
mkdir -p "$folder_temp_result/$sample_id/05-bakta"
mkdir -p "$folder_temp_result/$sample_id/06-quast"
mkdir -p "$folder_temp_result/$sample_id/07-busco"
mkdir -p "$folder_temp_result/$sample_id/08-mash"
mkdir -p "$folder_temp_result/$sample_id/09-sourmash"


# Specify log directory
log01="$folder_temp_result/$sample_id/logs/01-filter.log"
log02="$folder_temp_result/$sample_id/logs/02-miniasm.log"
log03="$folder_temp_result/$sample_id/logs/03-flye.log"
log04="$folder_temp_result/$sample_id/logs/04-medaka.log"
log05="$folder_temp_result/$sample_id/logs/05-bakta.log"
log06="$folder_temp_result/$sample_id/logs/06-quast.log"
log07="$folder_temp_result/$sample_id/logs/07-busco.log"
log08="$folder_temp_result/$sample_id/logs/08-mash.log"
log09="$folder_temp_result/$sample_id/logs/09-sourmash.log"

# 1. Filter worst reads via filtlong
raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
filtered_reads="$folder_temp_result/$sample_id/01-filtlong/01-filtered_reads.fastq.gz"

zsh 31-filter_reads.sh $raw_reads $filtered_reads &> $log01

# 2. Draft genome via miniasm
downsample_reads="$folder_temp_result/$sample_id/02-miniasm/02-downsampled_reads.fastq"
assemble_paf="$folder_temp_result/$sample_id/02-miniasm/02-miniasm_pilot.paf.gz"
assemble_gfa="$folder_temp_result/$sample_id/02-miniasm/02-miniasm_pilot.gfa"
assemble_fq="$folder_temp_result/$sample_id//02-miniasm/02-miniasm_pilot.fasta"

zsh 32-miniasm.sh $filtered_reads $downsample_reads $assemble_paf $assemble_gfa $assemble_fq &> $log02

# 3. Assemble genome via flye
filtered_reads="$folder_temp_result/Chang_Q5C_1/01-filtlong/01-filtered_reads.fastq.gz"
assemble_fq="$folder_temp_result/Chang_Q5C_1/02-miniasm/02-miniasm_pilot.fasta"
downsampled_reads2="$folder_temp_result/Chang_Q5C_1/03-flye/03-downsampled2_reads.fastq"
genome_size=$(grep -v ">" "$assemble_fq" | grep -E -o "G|C|T|A|N" | wc -l)
downsampled_reads2="$folder_temp_result/Chang_Q5C_1/03-flye/03-downsampled2_reads.fastq"
assembly_folder="$folder_temp_result/Chang_Q5C_1/03-flye/"

zsh 33-flye.sh $filtered_reads $assemble_fq $downsampled_reads2 $genome_size $downsampled_reads2 $assembly_folder | tee $log03

# 4. Polish genome via medaka
raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
draft_assembly="$folder_temp_result/Chang_Q5C_1/03-flye/assembly.fasta"
medaka_folder="$folder_temp_result/Chang_Q5C_1/04-medaka"

zsh 34-medaka.sh $raw_reads $draft_assembly $medaka_folder | tee $log04

# 5. Annotate genome via bakta
bakta_database="/Users/cychang/bioinformatics/bakta/db" # This database is mandatory and must be downloaded before annotation
medaka_consensus="$folder_temp_result/$sample_id/04-medaka/consensus.fasta"
bakta_folder="$folder_temp_result/$sample_id/05-bakta"

zsh 35-bakta.sh $bakta_database $medaka_consensus $bakta_folder | tee $log05

# 6. Check genome quality via quast
medaka_consensus="$folder_temp_result/$sample_id/04-medaka/consensus.fasta"
quast_folder="$folder_temp_result/$sample_id/06-quast"

zsh 36-quast.sh $medaka_consensus $quast_folder | tee $log06

# 7. Check genome quality via busco
medaka_consensus="$folder_temp_result/$sample_id/04-medaka/consensus.fasta"
busco_folder="$folder_temp_result/$sample_id/07-busco"

zsh 37-busco.sh $medaka_consensus $busco_folder | tee $log07

# 8. Identify taxonomy via mash
medaka_consensus="$folder_temp_result/$sample_id/04-medaka/consensus.fasta"
mash_sketch="$folder_temp_result/$sample_id/08-mash/consensus.fasta"
refseq_db="/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh"
mash_distance="$folder_temp_result/$sample_id/08-mash/distances.tab"

zsh 38-mash.sh $medaka_consensus $mash_sketch $refseq_db $mash_distance | tee $log08

# 8. Identify taxonomy via sourmash
medaka_consensus="$folder_temp_result/$sample_id/04-medaka/consensus.fasta"
sourmash_dir="$folder_temp_result/$sample_id/09-sourmash"
genbank_db="/Users/cychang/bioinformatics/sourmash/genbank-k31.lca.json.gz"
sourmash_sig="$folder_temp_result/$sample_id/09-sourmash/consensus.fasta.sig"

zsh 39-sourmash.sh $medaka_consensus $sourmash_dir $genbank_db $sourmash_sig | tee $log09


# Remove all temporary files











