# This rule runs nanoplot to qc nanopore raw reads
#genome_ids = ["g2", "g3"]

print(genome_ids)

# rule all:
#     input:
#         expand(
#             f"{config['folder_genomics']}/test/out/{{gi}}.txt",
#             gi = genome_ids, 
#         )

wildcards = glob_wildcards(f"{config['folder_genomics']}/test/{{sample}}.txt")
print(wildcards.sample)
SAMPLES = wildcards.sample

rule all:
    input:
        expand("results/{sample}_processed.txt", sample=SAMPLES)

rule test:
    input:
         "data/{sample}_data.txt"
    output:
        "results/{sample}_processed.txt"
    shell:
        "echo {wildcards.sample} > {output}"


    # echo $genome_ids[$i]
    # raw_reads=$folder_genomics/raw_reads/$genome_ids[$i].fastq.gz
    # filtered_reads=$folder_genomics/assembly/$genome_ids[$i]/filtered_reads.fastq.gz
    # mkdir -p $folder_genomics/assembly/$genome_ids[$i]
    # 
    # # Remove the worst 5% reads via filtlong
    # zsh 01a-filter_reads.sh \
    #     $raw_reads \
    #     $filtered_reads
