#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastqc } from './modules/fastqc.nf'
include { multiqc } from './modules/multiqc.nf'

workflow {
  Channel.fromFilePairs( "${params.run_dir}/Data/Intensities/BaseCalls/*_{R1,R2}_*.fastq.gz" ).set{ ch_fastq }

  main:
    fastqc(ch_fastq)
    multiqc(fastqc.out.map{ it -> [it[1], it[2]] }.collect())
}