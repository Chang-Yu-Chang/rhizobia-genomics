#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.raw_reads = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/raw_reads/Chang_Q5C_1.fastq.gz"
params.outdir = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics"
params.workdir = '/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/'

// Set custom work directory
// config. { workDir = params.workdir }

// Filter the worst 5% reads via filtlong
process filter_reads {
    input:
    path raw_reads

    output:
    path 'filtered_reads.fastq.gz'

    script:
    """
    #!/usr/bin/env zsh
    source ~/.zshrc

    # This throws away the worst 5% reads
    conda activate
    mamba activate filtlong

    filtlong --keep_percent 95 $raw_reads | gzip > filtered_reads.fastq.gz
    # `--keep_percent 95` throw out the worst 5% of reads

    """
}

// This extract the filtered reads into a txt file
process extract_reads {
    input:
    path filtered_reads

    publishDir '/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/filtered_reads'

    output:
    path 'filtered_reads.txt'

    script:
    """
    #!/usr/bin/env zsh
    source ~/.zshrc

    # This extracts the filtered reads into a txt file
    mamba activate bioawk

    bioawk -c fastx '{print \$name \$qual, length(\$seq)}' $filtered_reads > filtered_reads.txt
    """
}


workflow {
    // def proteins = Channel.fromPath( '/some/path/*.fa' )
    // def query_ch = Channel.fromPath(params.raw_reads)
    // def query_ch = Channel.fromPath("/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/raw_reads/*.fastq.gz")
    // filter_reads(params.raw_reads) | view
    // filter_reads(query_ch) | extract_reads
}



















