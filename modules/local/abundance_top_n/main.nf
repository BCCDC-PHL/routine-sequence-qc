process abundance_top_n {

    tag { meta.id + " / " + taxonomic_level }

    cpus 1
    executor 'local'
    errorStrategy 'ignore'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(_), path(bracken_abundances), val(taxonomic_level)

    output:
    tuple val(meta), path("${meta.id}_${taxonomic_level}_top_*.tsv"), val(taxonomic_level)

    script:
    top_n = taxonomic_level == 'Genus' ? '3' : '5'
    taxonomic_level_char = taxonomic_level.substring(0,1)
    """
    bracken_top_n_linelist.py ${bracken_abundances} -n ${top_n} -s ${meta.id} -l ${taxonomic_level_char} > ${meta.id}_${taxonomic_level}_top_${top_n}.tsv
    """
}
