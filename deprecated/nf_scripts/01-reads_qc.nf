#!/usr/bin/env nextflow

// Script parameters
params.raw_reads = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/raw_reads/*.fastq.gz"

// Filter the worst 5% reads via filtlong
process filter_reads {
    conda 'conda_envs/filtlong.yml'

    input:
    path raw_reads

    output:
    path 'filtered_reads'

    script:
    """
    filtlong --keep_percent 95 $raw_reads | gzip > filtered_reads
    """
    // --keep_percent 95 throw out the worst 5% of reads
}

// This script performs quality control on the raw nanapore whole-genome long reads
process extract_reads {
    conda 'conda_envs/bioawk.yml'

    input:
    path filtered_reads

    output:
    path "filtered_reads_*.txt"
    script:
    """
    bioawk -c fastx '{print \$name, \$qual, length(\$seq)}' $filtered_reads > filtered_reads_*.txt
    """
}



workflow {
    rr_ch = Channel.fromPath(params.raw_reads)
    filtered_reads = filter_reads(rr_ch)
    filtered_reads.view()

    extracted_reads = extract_reads(filtered_reads)
    extracted_reads.view()
    // def proteins = Channel.fromPath( '/some/path/*.fa' )
    // def query_ch = Channel.fromPath(params.raw_reads)
    // def query_ch = Channel.fromPath("/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/raw_reads/*.fastq.gz")
    // filter_reads(params.raw_reads) | view
    // filter_reads(query_ch) | extract_reads
}



















