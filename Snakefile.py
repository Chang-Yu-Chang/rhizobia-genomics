# Import Snakemake module
rule all:
    input:
        expand("/path/to/this/output/assembly/{sample}_assembly.fasta", sample=genomes)

# List of genome names
genomes = ["genome1", "genome2", "genome3", "genome4", "genome5",
           "genome6", "genome7", "genome8", "genome9", "genome10",
           "genome11", "genome12", "genome13", "genome14", "genome15",
           "genome16", "genome17", "genome18", "genome19", "genome20"]

# Mamba environments for each tool
conda: "envs/filtlong.yaml"
rule filtlong:
    input:
        "/path/to/this/input/{sample}.fasta"
    output:
        "/path/to/this/output/filtered/{sample}_filtered.fasta"
    shell:
        "filtlong --input {input} --output {output}"

conda: "envs/miniasm.yaml"
rule miniasm:
    input:
        "/path/to/this/output/filtered/{sample}_filtered.fasta"
    output:
        "/path/to/this/output/miniasm/{sample}_assembly.gfa"
    shell:
        "miniasm -f {input} > {output}"

conda: "envs/flye.yaml"
rule flye:
    input:
        "/path/to/this/output/miniasm/{sample}_assembly.gfa"
    output:
        "/path/to/this/output/flye/{sample}_assembly.fasta"
    shell:
        "flye --assembly {input} --out-dir /path/to/this/output/flye/{sample} --threads 4"

conda: "envs/medaka.yaml"
rule medaka:
    input:
        "/path/to/this/output/flye/{sample}_assembly.fasta"
    output:
        "/path/to/this/output/medaka/{sample}_polished.fasta"
    shell:
        "medaka_consensus -i {input} -d {input} -o {output}"

conda: "envs/bakta.yaml"
rule bakta:
    input:
        "/path/to/this/output/medaka/{sample}_polished.fasta"
    output:
        "/path/to/this/output/bakta/{sample}_annotation.gff"
    shell:
        "bakta {input} -o {output}"

conda: "envs/quast.yaml"
rule quast:
    input:
        "/path/to/this/output/medaka/{sample}_polished.fasta"
    output:
        "/path/to/this/output/quast/{sample}_quast_report"
    shell:
        "quast {input} -o /path/to/this/output/quast -r reference.fasta"

conda: "envs/busco.yaml"
rule busco:
    input:
        "/path/to/this/output/medaka/{sample}_polished.fasta"
    output:
        "/path/to/this/output/busco/{sample}_busco_results"
    shell:
        "busco -i {input} -o /path/to/this/output/busco/{sample} -l busco_dataset -m geno"

conda: "envs/mash.yaml"
rule mash:
    input:
        "/path/to/this/output/medaka/{sample}_polished.fasta"
    output:
        "/path/to/this/output/mash/{sample}.msh"
    shell:
        "mash sketch -o /path/to/this/output/mash/{sample} {input}"

conda: "envs/sourmash.yaml"
rule sourmash:
    input:
        "/path/to/this/output/mash/{sample}.msh"
    output:
        "/path/to/this/output/sourmash/{sample}_sourmash_results.txt"
    shell:
        "sourmash compare -k 21 -o /path/to/this/output/sourmash/{sample} {input}"

# Define a sample list
samples = ["sample1", "sample2", "sample3"]

# Create a list of input files for each sample
rule expand_samples:
    output:
        expand("/path/to/this/output/sourmash/{sample}_sourmash_results.txt", sample=samples)
