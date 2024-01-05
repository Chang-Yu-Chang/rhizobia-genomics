#!/usr/bin/env nextflow



// Filter the worst 5% reads via filtlong
process filter_worst {
    input:
    path raw_reads

    output:
    path 'filterd_reads'

    script:
    """
    zsh ../shell_scripts/01a-filter_reads.sh \
        $raw_reads \
        $filtered_reads
    """
}


workflow {
    filter_worst(params.raw_reads)
}
