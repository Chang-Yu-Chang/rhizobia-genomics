#!/usr/bin/env nextflow

//params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
//params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
//params.multiqc = "$projectDir/multiqc"
params {
  inputdir = 'path/to/your/input/'
  outdir = "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics"
  workdir = 'path/to/custom_work_dir'
}

// Set custom work directory
nextflow.config = [
  'workDir': params.workdir
]

log.info """\

    Bacterial Genome Assembly - NF Pipeline
    =======================================
    raw_reads    : ${params.raw_reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)
