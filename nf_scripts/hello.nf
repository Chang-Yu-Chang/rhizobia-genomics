#!/usr/bin/env nextflow

params.greeting = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

process splitLetters {

    input:
    val x

    output:
    path 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}


process convertUpper {

    input:
    path y

    output:
    stdout
    val y

    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

workflow my_pipeline2 {
    take:
        greeting_ch
    main:
        splitLetters(greeting_ch)
        convertUpper( splitLetters.out.flatten() )
        convertUpper.out[0].view { it }
}
