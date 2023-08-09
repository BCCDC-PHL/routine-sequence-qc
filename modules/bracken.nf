process bracken {

    tag { sample_id + " / " + taxonomic_level }

    errorStrategy 'ignore'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_${taxonomic_level}_bracken*", mode: 'copy'

    cpus 2

    input:
      tuple val(sample_id), path(kraken2_report), path(bracken_db), path(sample_sheet_json), val(taxonomic_level)

    output:
      tuple val(sample_id), path("${sample_id}_${taxonomic_level}_bracken.txt"), path("${sample_id}_${taxonomic_level}_bracken_abundances.tsv"), val(taxonomic_level), emit: unadjusted
      tuple val(sample_id), path("${sample_id}_${taxonomic_level}_bracken_adjusted.txt"), path("${sample_id}_${taxonomic_level}_bracken_abundances_adjusted.tsv"), val(taxonomic_level), emit: adjusted

    script:
    taxonomic_level_char = taxonomic_level.substring(0,1)
    """
    bracken -d ${bracken_db} \
      -i ${kraken2_report} \
      -w ${sample_id}_${taxonomic_level}_bracken.txt \
      -o ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted.tsv \
      -r \$(get_read_length.py ${sample_sheet_json}) \
      -l ${taxonomic_level_char}

    head -n 1 ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted.tsv > bracken_abundances_header.tsv
    tail -n+2 ${sample_id}_${taxonomic_level}_bracken_abundances_unsorted.tsv | sort -t \$'\\t' -nrk 7,7 > ${sample_id}_${taxonomic_level}_bracken_abundances_data.tsv
    cat bracken_abundances_header.tsv ${sample_id}_${taxonomic_level}_bracken_abundances_data.tsv > ${sample_id}_${taxonomic_level}_bracken_abundances.tsv

    adjust_for_unclassified_reads.py \
      --kraken-report ${kraken2_report} \
      --bracken-report ${sample_id}_${taxonomic_level}_bracken.txt \
      --bracken-abundances ${sample_id}_${taxonomic_level}_bracken_abundances.tsv \
      --adjusted-report ${sample_id}_${taxonomic_level}_bracken_adjusted.txt \
      --adjusted-abundances ${sample_id}_${taxonomic_level}_bracken_abundances_adjusted.tsv
    """
}


process abundance_top_n {

    tag { sample_id + " / " + taxonomic_level }

    errorStrategy 'ignore'

    executor 'local'

    cpus 1

    input:
      tuple val(sample_id), path(_), path(bracken_abundances), val(taxonomic_level)

    output:
      tuple val(sample_id), path("${sample_id}_${taxonomic_level}_top_*.tsv"), val(taxonomic_level)

    script:
    def top_n = taxonomic_level == 'Genus' ? '3' : '5'
    """
    bracken_top_n_linelist.py ${bracken_abundances} -n ${top_n} -s ${sample_id} > ${sample_id}_${taxonomic_level}_top_${top_n}.tsv
    """
}
