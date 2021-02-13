process bracken {

    tag { sample_id + " / " + taxonomic_level }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_${taxonomic_level}_bracken*", mode: 'copy'

    cpus 2

    input:
      tuple val(sample_id), path(kraken2_report), path(bracken_db), val(taxonomic_level)

    output:
      tuple val(sample_id), path("${sample_id}_${taxonomic_level}_bracken.txt"), path("${sample_id}_${taxonomic_level}_bracken_abundances.tsv"), val(taxonomic_level)

    script:
    taxonomic_level_char = taxonomic_level.substring(0,1)
    // grep/awk command adds a single space to the front of each line, for compatibility with MultiQC kraken module
    // MultiQC uses the following regex on the first two lines of a file to identify it as a kraken output:
    // '^\s{1,2}(\d{1,2}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKPCOFGS-]{1,3})\t(\d+)\s+(.+)'
    """
    bracken -d ${bracken_db} -i ${kraken2_report} -w ${sample_id}_${taxonomic_level}_bracken_tmp.txt -o ${sample_id}_${taxonomic_level}_bracken_abundances.tsv -r 250 -l ${taxonomic_level_char}
    sed 's/100\\.00/99\\.99/' ${sample_id}_${taxonomic_level}_bracken_tmp.txt | awk 'NR != 2' | awk '{print " ", \$0}' > ${sample_id}_${taxonomic_level}_bracken_tmp2.txt
    echo -e "  0.01\\t1\\t1\\tU\\t1\\tunclassified" > unclassified_placeholder.tsv
    cat unclassified_placeholder.tsv ${sample_id}_${taxonomic_level}_bracken_tmp2.txt > ${sample_id}_${taxonomic_level}_bracken.txt
    """
}

process abundance_top_n {

    tag { sample_id + " / " + taxonomic_level }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_${taxonomic_level}_top_*.tsv", mode: 'copy'

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
