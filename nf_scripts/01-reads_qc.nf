#!/usr/bin/env nextflow

// Script parameters
params.raw_reads = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/raw_reads/Chang_Q5C_1.fastq.gz"

// Filter the worst 5% reads via filtlong
process filter_reads {
    container 'filtlong'

    input:
    path raw_reads

    output:
    path 'filtered_reads'

    script:
    """
    mamba init 
    mamba activate filtlong
    filtlong --keep_percent 95 $raw_reads | gzip > filtered_reads
    # `--keep_percent 95` throw out the worst 5% of reads

    """
}

workflow {
    filtered_reads = filter_reads(params.raw_reads)
    filtered_reads.view()
    // def proteins = Channel.fromPath( '/some/path/*.fa' )
    // def query_ch = Channel.fromPath(params.raw_reads)
    // def query_ch = Channel.fromPath("/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/raw_reads/*.fastq.gz")
    // filter_reads(params.raw_reads) | view
    // filter_reads(query_ch) | extract_reads
}



















