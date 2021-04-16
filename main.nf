#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastqc } from './modules/fastqc.nf'
include { multiqc } from './modules/multiqc.nf'
include { interop_summary } from './modules/interop.nf'
include { parse_sample_sheet } from './modules/sample-sheet.nf'
include { kraken2 } from './modules/kraken2.nf'
include { bracken } from './modules/bracken.nf'
include { abundance_top_n } from './modules/bracken.nf'

workflow {
  ch_fastq = Channel.fromFilePairs( "${params.run_dir}/Data/Intensities/BaseCalls/*_{R1,R2}_*.fastq.gz" )
  ch_sample_sheet = Channel.fromPath( "${params.run_dir}/SampleSheet*.csv" )
  ch_multiqc_config = Channel.fromPath( "${projectDir}/assets/multiqc_config_base.yaml" )
  ch_run_dir = Channel.fromPath(params.run_dir)
  ch_run_id = Channel.fromPath(params.run_dir).map{ it -> it.baseName }
  ch_kraken2_db = Channel.fromPath(params.kraken2_db)
  ch_bracken_db = Channel.fromPath(params.bracken_db)
  ch_taxonomic_levels = Channel.from('Genus', 'Species')

  main:
    
    interop_summary(ch_run_id.combine(ch_run_dir))

    ch_sample_sheet.view()

    parse_sample_sheet(ch_run_id.combine(ch_sample_sheet))

    fastqc(ch_fastq)

    kraken2(ch_fastq.combine(ch_kraken2_db))
    bracken(kraken2.out.combine(ch_bracken_db).combine(parse_sample_sheet.out).combine(ch_taxonomic_levels))

    abundance_top_n(bracken.out)

    abundance_top_n.out.filter{ it[2] == 'Genus' }.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "top_3_abundances_genus.csv", storeDir: "${params.outdir}/abundance_top_n")
    abundance_top_n.out.filter{ it[2] == 'Species' }.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "top_5_abundances_species.csv", storeDir: "${params.outdir}/abundance_top_n")
    
    ch_fastqc_collected = fastqc.out.map{ it -> [it[1], it[2]] }.collect()
    ch_bracken_species_multiqc_collected = bracken.out.filter{ it[4] == 'Species' }.map{ it -> it[2] }.collect()
    ch_all_qc_outputs = interop_summary.out.map{ it -> it.drop(1) }.combine(ch_fastqc_collected).combine(ch_bracken_species_multiqc_collected)
    
    ch_all_qc_outputs_with_run_id = ch_run_id.combine(ch_all_qc_outputs).map{ it -> [it[0], it.drop(1)] }
    multiqc(ch_multiqc_config.combine(ch_all_qc_outputs_with_run_id))
}