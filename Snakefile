# Define a base directory path
folder_data = "/user/cychang/Dropbox/lab/local-adaptation/data/"

# # Import Snakemake module
rule all:
    input:
    expand("raw/assembly/{sample}_assembly.fasta", sample=genomes)

# List of genome names
#sample_batch_id = [i for i in ["Chang_Q5C_results"] for _ in range(10)] + ["Chang_Q5C_results_repeated"] + [i for i in #["Chang_Q5C_results"] for _ in range(6)] + ["Chang_Q5C_results_repeated", "Chang_Q5C_results"]
#sample_id = ["Chang_Q5C_" + str(i) for i in range(1,20)]

sample_batch_id = "Chang_Q5C_results"
sample_id = "Chang_Q5C_1"

# Create directory
rule create_directory:
    input:
        # Any dependencies or input files
    output:
        directory("{folder_data}/temp/plasmidsaurus/")
    threads: 8
    shell:
        """
        mkdir -p "{sample_batch_id}/test"
        """

# Mamba environments for each tool
#conda: "envs/filtlong.yaml"
rule filtlong:
    input:
        f"{folder_data}/raw/{{sample_batch_id}}/{{sample_id}}/reads/raw_reads.fasta.gz"
    output:
        f"{folder_data}/temp/plasmidsaurus/{{sample_batch_id}}/{{sample_id}}/01-filtlong/01-filtered_reads.fastq"
    shell:
        """
        cd
        source ~/.zshrc
        # Throw out the worst 5% reads
        conda activate
        mamba activate filtlong
        mamba env list
        filtlong --keep_percent 95 $raw_reads | gzip > $filtered_reads
        """
