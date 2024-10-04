# This rule runs nanoplot to qc nanopore raw reads
genome_ids = ["g2", "g3"]

rule nanoplot:
    input:
        lambda wildcards: f"{config['folder_genomics']}/raw_reads/{wildcards.gi}.fastq.gz",
    output:
        f"{config['folder_genomics']}/nanoplot/{{gi}}/report.html",
    log:
        f"{config['folder_genomics']}/nanoplot/{{gi}}/log",
    conda:
        "../envs/nanoplot.yaml"
    shell:
        """
        NanoPlot -t 20 --fastq {input} --loglength -o {output} --plots dot --verbose >> {log} 2>&1
        """


    # echo $genome_ids[$i]
    # raw_reads=$folder_genomics/raw_reads/$genome_ids[$i].fastq.gz
    # filtered_reads=$folder_genomics/assembly/$genome_ids[$i]/filtered_reads.fastq.gz
    # mkdir -p $folder_genomics/assembly/$genome_ids[$i]
    # 
    # # Remove the worst 5% reads via filtlong
    # zsh 01a-filter_reads.sh \
    #     $raw_reads \
    #     $filtered_reads
