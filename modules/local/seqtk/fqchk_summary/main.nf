process seqtk_fqchk_summary {

    tag { meta.id }

    cpus 1
    conda "${moduleDir}/environment.yml"

    executor 'local'

    input:
    tuple val(meta), path(seqtk_fqchk_output_all_positions), path(seqtk_fqchk_output_by_position)

    output:
    tuple val(meta), path("${meta.id}_seqtk_fqchk_summary.csv")

    script:
    """
    summarize_seqtk_fqchk.py ${seqtk_fqchk_output_all_positions} --sample-id ${meta.id} > ${meta.id}_seqtk_fqchk_summary.csv
    """
}
