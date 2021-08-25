#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastqc } from './modules/fastqc.nf'
include { seqtk_fqchk } from './modules/seqtk.nf'
include { seqtk_fqchk_summary } from './modules/seqtk.nf'
include { multiqc } from './modules/multiqc.nf'
include { interop_summary } from './modules/interop.nf'
include { parse_sample_sheet } from './modules/sample-sheet.nf'
include { kraken2 } from './modules/kraken2.nf'
include { bracken } from './modules/bracken.nf'
include { abundance_top_n } from './modules/bracken.nf'
include { mash_sketch } from './modules/mash.nf'
include { mash_sketch_summary } from './modules/mash.nf'
include { combine_qc_stats } from './modules/combine_qc_stats.nf'
include { copy_to_nginx } from './modules/copy_to_nginx.nf'

if (params.instrument_type == "miseq") {
  fastq_subdir = "Data/Intensities/BaseCalls"
} else if (params.instrument_type == "nextseq") {
  fastq_subdir = "Analysis/1/Data/fastq"
} else {
  println("Unsupported instrument type: ${params.instrument_type}")
  System.exit(1)
}

workflow {
  ch_fastq = Channel.fromFilePairs( "${params.run_dir}/${fastq_subdir}/*_{R1,R2}*.fastq.gz" )
  ch_sample_sheet = Channel.fromPath( "${params.run_dir}/SampleSheet*.csv" )
  ch_multiqc_config = Channel.fromPath( "${projectDir}/assets/multiqc_config_base.yaml" )
  ch_run_dir = Channel.fromPath(params.run_dir)
  ch_run_id = Channel.fromPath(params.run_dir).map{ it -> it.baseName }
  ch_kraken2_db = Channel.fromPath(params.kraken2_db)
  ch_bracken_db = Channel.fromPath(params.bracken_db)
  ch_taxonomic_levels = Channel.from('Genus', 'Species')

  // where the nginx server directory lives
  ch_nginx_server_location = Channel.fromPath("/srv/www/nginx/html2/csv-to-html-table/")

  main:
    
    interop_summary(ch_run_id.combine(ch_run_dir))

    parse_sample_sheet(ch_run_id.combine(ch_sample_sheet))

    fastqc(ch_fastq)

    seqtk_fqchk(ch_fastq)
    seqtk_fqchk_summary(seqtk_fqchk.out).map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "seqtk_fqchk_summary.csv", storeDir: "${params.outdir}/seqtk_fqchk_summary")

    mash_sketch(ch_fastq)
    mash_sketch_summary(mash_sketch.out).map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "mash_sketch_summary.csv", storeDir: "${params.outdir}/mash_sketch_summary")

    kraken2(ch_fastq.combine(ch_kraken2_db))
    bracken(kraken2.out.combine(ch_bracken_db).combine(parse_sample_sheet.out).combine(ch_taxonomic_levels))

    abundance_top_n(bracken.out)

    abundance_top_n.out.filter{ it[2] == 'Genus' }.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "top_3_abundances_genus.csv", storeDir: "${params.outdir}/abundance_top_n")
    abundance_top_5 = abundance_top_n.out.filter{ it[2] == 'Species' }.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "top_5_abundances_species.csv", storeDir: "${params.outdir}/abundance_top_n")

    ch_fastqc_collected = fastqc.out.map{ it -> [it[1], it[2]] }.collect()
    ch_bracken_species_multiqc_collected = bracken.out.filter{ it[4] == 'Species' }.map{ it -> it[2] }.collect()
    ch_all_qc_outputs = interop_summary.out.map{ it -> it.drop(1) }.combine(ch_fastqc_collected).combine(ch_bracken_species_multiqc_collected)

    basic_qc_stats = combine_qc_stats(abundance_top_n.out.filter{ it[2] == 'Species' }.map{ it -> [it[0], it[1]] }.join(seqtk_fqchk_summary.out.map{ it -> [it[0], it[1]] }).join(mash_sketch_summary.out)).map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "basic_qc_stats.csv", storeDir: "${params.outdir}/basic_qc_stats")

    ch_all_qc_outputs_with_run_id = ch_run_id.combine(ch_all_qc_outputs).map{ it -> [it[0], it.drop(1)] }
    multiqc_file = multiqc(ch_multiqc_config.combine(ch_all_qc_outputs_with_run_id))

    // not the best way to grab the html file, consider using regex for .html instead
    multiqc_file = multiqc_file.map{ it -> it[0]}
    copy_to_nginx(ch_run_id.combine(ch_nginx_server_location).combine(multiqc_file).combine(basic_qc_stats).combine(abundance_top_5))
}
