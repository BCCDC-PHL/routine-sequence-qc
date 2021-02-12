process estimate_abundance {

    tag { sample_id + " / " + taxonomic_level }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_*racken*.txt", mode: 'copy'

    cpus 8

    input:
      tuple val(sample_id), path(kraken2_report), path(bracken_db), val(taxonomic_level)

    output:
      tuple val(sample_id), path("${sample_id}_${taxonomic_level}_kraken2_updated.txt"), path("${sample_id}_${taxonomic_level}_bracken_report.txt")

    script:

    """
    est_abundance.py  -i ${kraken2_report} -k ${bracken_db}/database150mers.kmer_distrib -o ${sample_id}_${taxonomic_level}_kraken2_updated.txt -l ${taxonomic_level.substring(0,1)} --out-report ${sample_id}_${taxonomic_level}_bracken_report.txt
    """
}
