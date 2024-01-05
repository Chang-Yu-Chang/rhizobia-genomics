#!/usr/bin/env nextflow


// test_run.nf

// Define a process named 'appendTimestamp'
process appendTimestamp {

    // Define input file parameter
    input:
    file input_file

    // Define output file parameter
    output:
    file "output_${input_file.baseName}_${timestamp()}.txt"

    // Define the script to execute
    script:
    """
    #!/bin/bash
    date > ${output}
    echo "Content of input file:"
    cat ${input_file} >> ${output}
    """
}

// Define the workflow
workflow {

    // Specify an input file for the 'appendTimestamp' process
    inputFile = file("input.txt")

    // Execute the 'appendTimestamp' process with the specified input file
    appendTimestamp(input_file: inputFile)
}
