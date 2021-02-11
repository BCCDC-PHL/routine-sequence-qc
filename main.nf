#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastqc } from './modules/fastqc.nf'
include { multiqc } from './modules/multiqc.nf'
include { interop_summary } from './modules/interop.nf'

workflow {
  ch_fastq = Channel.fromFilePairs( "${params.run_dir}/Data/Intensities/BaseCalls/*_{R1,R2}_*.fastq.gz" )
  ch_run_dir = Channel.fromPath(params.run_dir)
  ch_run_id = Channel.fromPath(params.run_dir).map{ it -> it.baseName }

  main:
    interop_summary(ch_run_id.combine(ch_run_dir))

    fastqc(ch_fastq)
    
    // Line below is just composing the run_id and the list of fastqc_outdirs into a new list. There must be a better way(?)
    // "run_id" + ["fastqc_outdir1", "fastqc_outdir2", ...] => ["run_id", ["fastqc_outdir1", "fastqc_outdir2"]]
    multiqc(ch_run_id.combine(fastqc.out.map{ it -> [it[1], it[2]] }.collect()).map{ it -> [it[0], it.drop(1)] })

}