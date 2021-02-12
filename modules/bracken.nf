process estimate_abundance {

    tag { sample_id + " / " + taxonomic_level }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_${taxonomic_level}_bracken.tsv", mode: 'copy'

    cpus 2

    input:
      tuple val(sample_id), path(kraken2_report), path(bracken_db), val(taxonomic_level)

    output:
      tuple val(sample_id), path("${sample_id}_${taxonomic_level}_bracken.tsv"), val(taxonomic_level)

    script:

    """
    est_abundance.py  -i ${kraken2_report} -k ${bracken_db}/database150mers.kmer_distrib -o ${sample_id}_${taxonomic_level}_bracken.tsv -l ${taxonomic_level.substring(0,1)} --out-report ${sample_id}_${taxonomic_level}_bracken_report.txt
    """
}

process abundance_top_n {

    tag { sample_id + " / " + taxonomic_level }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_${taxonomic_level}_top_*.tsv", mode: 'copy'

    executor 'local'

    cpus 1

    input:
      tuple val(sample_id), path(bracken_report), val(taxonomic_level)

    output:
      tuple val(sample_id), path("${sample_id}_${taxonomic_level}_top_*.tsv")

    script:
    def top_n = taxonomic_level == 'G' ? '3' : '5'
    """
    bracken_top_n_linelist.py ${bracken_report} -n ${top_n} -s ${sample_id} > ${sample_id}_${taxonomic_level}_top_${top_n}.tsv
    """
}
