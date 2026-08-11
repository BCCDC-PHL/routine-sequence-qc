process bracken {

    tag { meta.id + " / " + taxonomic_level }

    cpus 2
    errorStrategy 'ignore'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/bracken", pattern: "${meta.id}_${taxonomic_level}_bracken*", mode: 'copy'

    input:
    tuple val(meta), path(kraken2_report), path(bracken_db), path(sample_sheet_json), val(taxonomic_level)

    output:
    tuple val(meta), path("${meta.id}_${taxonomic_level}_bracken.txt"), path("${meta.id}_${taxonomic_level}_bracken_abundances.tsv"), val(taxonomic_level), emit: unadjusted
    tuple val(meta), path("${meta.id}_${taxonomic_level}_bracken_adjusted.txt"), path("${meta.id}_${taxonomic_level}_bracken_abundances_adjusted.tsv"), val(taxonomic_level), emit: adjusted

    script:
    taxonomic_level_char = taxonomic_level.substring(0,1)
    """
    bracken -d ${bracken_db} \
      -i ${kraken2_report} \
      -w ${meta.id}_${taxonomic_level}_bracken.txt \
      -o ${meta.id}_${taxonomic_level}_bracken_abundances_unsorted.tsv \
      -r \$(get_read_length.py ${sample_sheet_json}) \
      -l ${taxonomic_level_char}

    head -n 1 ${meta.id}_${taxonomic_level}_bracken_abundances_unsorted.tsv > bracken_abundances_header.tsv
    tail -n+2 ${meta.id}_${taxonomic_level}_bracken_abundances_unsorted.tsv | sort -t \$'\\t' -nrk 7,7 > ${meta.id}_${taxonomic_level}_bracken_abundances_data.tsv
    cat bracken_abundances_header.tsv ${meta.id}_${taxonomic_level}_bracken_abundances_data.tsv > ${meta.id}_${taxonomic_level}_bracken_abundances.tsv

    adjust_for_unclassified_reads.py \
      --kraken-report ${kraken2_report} \
      --bracken-report ${meta.id}_${taxonomic_level}_bracken.txt \
      --bracken-abundances ${meta.id}_${taxonomic_level}_bracken_abundances.tsv \
      --adjusted-report ${meta.id}_${taxonomic_level}_bracken_adjusted.txt \
      --adjusted-abundances ${meta.id}_${taxonomic_level}_bracken_abundances_adjusted.tsv
    """
}
