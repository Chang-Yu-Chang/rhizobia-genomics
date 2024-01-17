#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.sample_list = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/sample_list.txt"
params.raw_data_folder = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/test0"
params.intermediate_folder1 = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/test1"
params.intermediate_folder2 = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/test2"

samples = file(params.sample_list)


//
process test1 {
    input:
    val sample_id from samples

    // publishDir '$params.genomics/test1/'
    output:
    path "${sample_id}_1.txt" into x

    script:
    """
    #!/usr/bin/env zsh
    cat ${raw_data_folder}/${sample_id}.txt > ${intermediate_folder1}/${sample_id}_1.txt
    """
}

//
process test2 {
    input:
    val sample_id from samples
    path x_file from x

    // publishDir '$genomics/test2'
    output:
    path path "${sample_id}_2.txt"

    script:
    """
    #!/usr/bin/env zsh
    cat  ${x_file} > ${intermediate_folder2}/${sample_id}_2.txt
    """
}


workflow {
    // Channel of sample IDs
    //sample_ch = Channel.from(samples)
     // Process each sample in parallel
    tes1(sample)
    tes2(sample)

    publishDir path: params.intermediate_folder1, pattern: '*_1.txt', mode: 'copy'
    publishDir path: params.intermediate_folder2, pattern: '*_2.txt', mode: 'copy'

}


